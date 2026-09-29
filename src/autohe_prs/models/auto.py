"""Optional boosted-tree models with a dependency-free ridge fallback."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from pathlib import Path
import pickle
from typing import Any

from .features import feature_names, feature_vector
from .train import (
    CostModelBundle,
    REGRESSION_TARGETS,
    _binary_label,
    _classification_metrics,
    _metrics,
    grouped_split,
    measured_training_rows,
    train_cost_models,
    usable_binary_label,
    usable_numeric,
)


@dataclass
class EstimatorRegressor:
    target: str
    estimator: Any
    log_target: bool = True
    model_feature_names: tuple[str, ...] = ()

    def predict(self, row: dict[str, Any]) -> float:
        value = float(
            self.estimator.predict(
                [feature_vector(row, self.model_feature_names or feature_names())]
            )[0]
        )
        return max(0.0, math.expm1(value)) if self.log_target else value


@dataclass
class EstimatorClassifier:
    estimator: Any
    model_feature_names: tuple[str, ...] = ()

    def predict_probability(self, row: dict[str, Any]) -> float:
        if hasattr(self.estimator, "predict_proba"):
            return float(
                self.estimator.predict_proba(
                    [feature_vector(row, self.model_feature_names or feature_names())]
                )[0][1]
            )
        return float(
            self.estimator.predict(
                [feature_vector(row, self.model_feature_names or feature_names())]
            )[0]
        )


@dataclass
class EstimatorCostModelBundle:
    regressors: dict[str, EstimatorRegressor]
    success_classifier: EstimatorClassifier | None
    error_feasibility_classifier: EstimatorClassifier | None
    metrics: dict[str, dict[str, float]]
    split: dict[str, Any]
    model_type: str
    target_scopes: dict[str, str] = field(default_factory=dict)

    def predict(self, row: dict[str, Any]) -> dict[str, float]:
        values = {target: model.predict(row) for target, model in self.regressors.items()}
        if self.success_classifier:
            values["execution_success_probability"] = (
                self.success_classifier.predict_probability(row)
            )
        if self.error_feasibility_classifier:
            values["error_constraint_feasibility_probability"] = (
                self.error_feasibility_classifier.predict_probability(row)
            )
        return values

    def save(self, output_dir: str | Path) -> None:
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        with (output / "cost_model.pkl").open("wb") as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)
        (output / "cost_model.json").write_text(
            json.dumps(
                {
                    "model_type": self.model_type,
                    "feature_names": feature_names(),
                    "targets": sorted(self.regressors),
                    "target_scopes": self.target_scopes,
                    "has_success_classifier": self.success_classifier is not None,
                    "has_error_feasibility_classifier": (
                        self.error_feasibility_classifier is not None
                    ),
                    "metrics": self.metrics,
                    "split": self.split,
                    "pickle_warning": "Load only models produced by a trusted AutoHE-PRS run.",
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )


def load_cost_model(path: str | Path) -> CostModelBundle | EstimatorCostModelBundle:
    path = Path(path)
    directory = path if path.is_dir() else path.parent
    pickle_path = directory / "cost_model.pkl"
    if pickle_path.exists():
        with pickle_path.open("rb") as handle:
            value = pickle.load(handle)
        if not isinstance(value, EstimatorCostModelBundle):
            raise ValueError("unrecognized AutoHE-PRS estimator bundle")
        return value
    return CostModelBundle.load(path)


def _factories(backend: str, seed: int) -> tuple[str, Any, Any]:
    requested = backend.lower()
    if requested in {"auto", "xgboost"}:
        try:
            from xgboost import XGBClassifier, XGBRegressor  # type: ignore[import-not-found]
        except ImportError:
            if requested == "xgboost":
                raise RuntimeError("xgboost backend requested but xgboost is not installed")
        else:
            return (
                "xgboost",
                lambda: XGBRegressor(
                    n_estimators=300, max_depth=6, learning_rate=0.05,
                    subsample=0.9, colsample_bytree=0.9, random_state=seed, n_jobs=1,
                ),
                lambda: XGBClassifier(
                    n_estimators=300, max_depth=6, learning_rate=0.05,
                    subsample=0.9, colsample_bytree=0.9, random_state=seed, n_jobs=1,
                ),
            )
    if requested in {"auto", "lightgbm"}:
        try:
            from lightgbm import LGBMClassifier, LGBMRegressor  # type: ignore[import-not-found]
        except ImportError:
            if requested == "lightgbm":
                raise RuntimeError("lightgbm backend requested but lightgbm is not installed")
        else:
            return (
                "lightgbm",
                lambda: LGBMRegressor(n_estimators=300, random_state=seed, n_jobs=1),
                lambda: LGBMClassifier(n_estimators=300, random_state=seed, n_jobs=1),
            )
    if requested in {"auto", "hist_gradient_boosting"}:
        try:
            from sklearn.ensemble import (  # type: ignore[import-not-found]
                HistGradientBoostingClassifier,
                HistGradientBoostingRegressor,
            )
        except ImportError:
            if requested == "hist_gradient_boosting":
                raise RuntimeError("scikit-learn backend requested but it is not installed")
        else:
            return (
                "hist_gradient_boosting",
                lambda: HistGradientBoostingRegressor(
                    max_iter=300, min_samples_leaf=5, random_state=seed
                ),
                lambda: HistGradientBoostingClassifier(
                    max_iter=300, min_samples_leaf=5, random_state=seed
                ),
            )
    if requested == "random_forest":
        try:
            from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
        except ImportError as exc:
            raise RuntimeError("random_forest requires scikit-learn") from exc
        return (
            "random_forest",
            lambda: RandomForestRegressor(n_estimators=400, random_state=seed, n_jobs=1),
            lambda: RandomForestClassifier(n_estimators=400, random_state=seed, n_jobs=1),
        )
    if requested == "mlp":
        try:
            from sklearn.neural_network import MLPClassifier, MLPRegressor
        except ImportError as exc:
            raise RuntimeError("mlp requires scikit-learn") from exc
        return (
            "mlp",
            lambda: MLPRegressor(
                hidden_layer_sizes=(64, 32), max_iter=1000, random_state=seed
            ),
            lambda: MLPClassifier(
                hidden_layer_sizes=(64, 32), max_iter=1000, random_state=seed
            ),
        )
    if requested not in {"auto", "ridge"}:
        raise ValueError(f"unsupported model backend: {backend}")
    raise ImportError


def train_auto_cost_models(
    rows: list[dict[str, Any]],
    *,
    backend: str = "auto",
    group_key: str = "workload_group",
    test_fraction: float = 0.2,
    seed: int = 20260728,
    held_out_groups: tuple[str, ...] | None = None,
) -> CostModelBundle | EstimatorCostModelBundle:
    rows = measured_training_rows(rows)
    if not rows:
        raise ValueError("cost-model training requires measured benchmark records")
    if backend.lower() == "ridge":
        return train_cost_models(
            rows, group_key=group_key, test_fraction=test_fraction, seed=seed,
            held_out_groups=held_out_groups,
        )
    try:
        model_type, regressor_factory, classifier_factory = _factories(backend, seed)
    except ImportError:
        return train_cost_models(
            rows, group_key=group_key, test_fraction=test_fraction, seed=seed,
            held_out_groups=held_out_groups,
        )

    successful = [row for row in rows if str(row.get("status")) == "ok"]
    train, test, test_groups = grouped_split(
        successful, group_key, test_fraction, seed, held_out_groups
    )
    if not train:
        raise ValueError("cost-model training requires successful training rows")
    regressors: dict[str, EstimatorRegressor] = {}
    metrics: dict[str, dict[str, float]] = {}
    target_scopes: dict[str, str] = {}
    for target in REGRESSION_TARGETS:
        target_train = [row for row in train if usable_numeric(row, target)]
        target_test = [row for row in test if usable_numeric(row, target)]
        if not target_train:
            continue
        if target == "key_bytes":
            scopes = {
                str(row.get("key_bytes_scope"))
                for row in target_train
                if row.get("key_bytes_scope")
            }
            if len(scopes) > 1:
                raise ValueError(
                    "key_bytes regression cannot mix incomparable scopes: "
                    + ", ".join(sorted(scopes))
                )
            if scopes:
                scope = next(iter(scopes))
                target_scopes[target] = scope
                target_test = [
                    row
                    for row in target_test
                    if str(row.get("key_bytes_scope")) == scope
                ]
        estimator = regressor_factory()
        estimator.fit(
            [feature_vector(row) for row in target_train],
            [math.log1p(max(0.0, float(row[target]))) for row in target_train],
        )
        model = EstimatorRegressor(target, estimator, True, feature_names())
        regressors[target] = model
        metrics[target] = _metrics(model, target_test or target_train, target)
    all_train, all_test, _ = grouped_split(
        rows, group_key, test_fraction, seed, held_out_groups
    )
    labels = [int(str(row.get("status")) == "ok") for row in all_train]
    classifier = None
    if len(set(labels)) > 1:
        estimator = classifier_factory()
        estimator.fit([feature_vector(row) for row in all_train], labels)
        classifier = EstimatorClassifier(estimator, feature_names())
    metrics["execution_success"] = _classification_metrics(
        classifier, all_test or all_train  # type: ignore[arg-type]
    )
    error_rows = [
        row for row in rows
        if usable_binary_label(row, "error_constraint_feasible")
    ]
    try:
        error_train, error_test, error_test_groups = grouped_split(
            error_rows, group_key, test_fraction, seed, held_out_groups
        )
    except ValueError:
        error_train, error_test, error_test_groups = error_rows, [], []
    error_labels = [
        int(_binary_label(row, "error_constraint_feasible"))
        for row in error_train
    ]
    error_classifier = None
    if len(set(error_labels)) > 1:
        estimator = classifier_factory()
        estimator.fit([feature_vector(row) for row in error_train], error_labels)
        error_classifier = EstimatorClassifier(estimator, feature_names())
    metrics["error_constraint_feasibility"] = _classification_metrics(
        error_classifier,  # type: ignore[arg-type]
        error_test or error_train,
        label_key="error_constraint_feasible",
    )
    return EstimatorCostModelBundle(
        regressors,
        classifier,
        error_classifier,
        metrics,
        {
            "group_key": group_key,
            "test_fraction": test_fraction,
            "seed": seed,
            "test_groups": test_groups,
            "train_rows": len(train),
            "test_rows": len(test),
            "error_feasibility_train_rows": len(error_train),
            "error_feasibility_test_rows": len(error_test),
            "error_feasibility_test_groups": error_test_groups,
        },
        model_type,
        target_scopes,
    )
