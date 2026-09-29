"""Dependency-free ridge/logistic cost models with group-held-out evaluation."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
import random
from typing import Any

from .features import feature_names, feature_vector


REGRESSION_TARGETS = (
    "evaluation_seconds",
    "total_seconds",
    "peak_memory_mb",
    "ciphertext_bytes",
    "key_bytes",
    "max_absolute_error",
    "max_relative_error",
)


def measured_training_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Exclude predictions/simulations from empirical target fitting."""
    return [
        row for row in rows
        if str(row.get("measurement_kind", "measured_unspecified")).startswith("measured")
        or row.get("measurement_kind") == "validation_only"
    ]


def usable_numeric(row: dict[str, Any], target: str) -> bool:
    try:
        return math.isfinite(float(row[target]))
    except (KeyError, TypeError, ValueError):
        return False


def usable_binary_label(row: dict[str, Any], label_key: str) -> bool:
    raw = row.get(label_key)
    if raw in {None, ""}:
        return False
    if isinstance(raw, float) and not math.isfinite(raw):
        return False
    return True


def _solve(matrix: list[list[float]], vector: list[float]) -> list[float]:
    n = len(vector)
    augmented = [row[:] + [value] for row, value in zip(matrix, vector)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        if abs(augmented[column][column]) < 1e-12:
            augmented[column][column] = 1e-12
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [augmented[row][-1] for row in range(n)]


@dataclass(frozen=True)
class LinearModel:
    target: str
    coefficients: tuple[float, ...]
    intercept: float
    log_target: bool
    model_feature_names: tuple[str, ...] = ()

    def predict(self, row: dict[str, Any]) -> float:
        raw = self.intercept + sum(
            weight * value
            for weight, value in zip(
                self.coefficients,
                feature_vector(row, self.model_feature_names or feature_names()),
            )
        )
        return max(0.0, math.expm1(raw)) if self.log_target else raw


@dataclass(frozen=True)
class LogisticModel:
    coefficients: tuple[float, ...]
    intercept: float
    model_feature_names: tuple[str, ...] = ()

    def predict_probability(self, row: dict[str, Any]) -> float:
        raw = self.intercept + sum(
            weight * value
            for weight, value in zip(
                self.coefficients,
                feature_vector(row, self.model_feature_names or feature_names()),
            )
        )
        if raw >= 0:
            return 1.0 / (1.0 + math.exp(-raw))
        exp_value = math.exp(raw)
        return exp_value / (1.0 + exp_value)


@dataclass(frozen=True)
class CostModelBundle:
    feature_names: tuple[str, ...]
    regressors: dict[str, LinearModel]
    success_classifier: LogisticModel | None
    error_feasibility_classifier: LogisticModel | None
    metrics: dict[str, dict[str, float]]
    split: dict[str, Any]
    model_type: str = "ridge_closed_form"
    target_scopes: dict[str, str] = field(default_factory=dict)

    def predict(self, row: dict[str, Any]) -> dict[str, float]:
        result = {target: model.predict(row) for target, model in self.regressors.items()}
        if self.success_classifier:
            result["execution_success_probability"] = self.success_classifier.predict_probability(row)
        if self.error_feasibility_classifier:
            result["error_constraint_feasibility_probability"] = (
                self.error_feasibility_classifier.predict_probability(row)
            )
        return result

    def save(self, output_dir: str | Path) -> None:
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        value = asdict(self)
        (output / "cost_model.json").write_text(
            json.dumps(value, indent=2, sort_keys=True), encoding="utf-8"
        )

    @staticmethod
    def load(path: str | Path) -> "CostModelBundle":
        path = Path(path)
        if path.is_dir():
            path = path / "cost_model.json"
        raw = json.loads(path.read_text(encoding="utf-8"))
        models = {
            target: LinearModel(
                target=value["target"],
                coefficients=tuple(value["coefficients"]),
                intercept=value["intercept"],
                log_target=value["log_target"],
                model_feature_names=tuple(
                    value.get("model_feature_names") or raw["feature_names"]
                ),
            )
            for target, value in raw["regressors"].items()
        }
        classifier_raw = raw.get("success_classifier")
        classifier = (
            LogisticModel(
                tuple(classifier_raw["coefficients"]),
                classifier_raw["intercept"],
                tuple(classifier_raw.get("model_feature_names") or raw["feature_names"]),
            )
            if classifier_raw else None
        )
        error_raw = raw.get("error_feasibility_classifier")
        error_classifier = (
            LogisticModel(
                tuple(error_raw["coefficients"]),
                error_raw["intercept"],
                tuple(error_raw.get("model_feature_names") or raw["feature_names"]),
            )
            if error_raw else None
        )
        target_scopes = dict(raw.get("target_scopes", {}))
        if "key_bytes" in models and "key_bytes" not in target_scopes:
            target_scopes["key_bytes"] = "legacy_context+public_key_proxy"
        return CostModelBundle(
            feature_names=tuple(raw["feature_names"]),
            regressors=models,
            success_classifier=classifier,
            error_feasibility_classifier=error_classifier,
            metrics=raw["metrics"],
            split=raw["split"],
            model_type=raw["model_type"],
            target_scopes=target_scopes,
        )


def _fit_ridge(rows: list[dict[str, Any]], target: str, alpha: float) -> LinearModel:
    x = [[1.0, *feature_vector(row)] for row in rows]
    y = [math.log1p(max(0.0, float(row[target]))) for row in rows]
    dimensions = len(x[0])
    gram = [[0.0] * dimensions for _ in range(dimensions)]
    rhs = [0.0] * dimensions
    for values, target_value in zip(x, y):
        for i in range(dimensions):
            rhs[i] += values[i] * target_value
            for j in range(dimensions):
                gram[i][j] += values[i] * values[j]
    for i in range(1, dimensions):
        gram[i][i] += alpha
    coefficients = _solve(gram, rhs)
    return LinearModel(
        target, tuple(coefficients[1:]), coefficients[0], True, feature_names()
    )


def _metrics(model: LinearModel, rows: list[dict[str, Any]], target: str) -> dict[str, float]:
    rows = [row for row in rows if usable_numeric(row, target)]
    if not rows:
        return {"count": 0.0}
    observed = [float(row[target]) for row in rows]
    predicted = [model.predict(row) for row in rows]
    errors = [prediction - actual for prediction, actual in zip(predicted, observed)]
    mae = sum(abs(value) for value in errors) / len(errors)
    rmse = math.sqrt(sum(value * value for value in errors) / len(errors))
    mape_values = [abs(error / actual) for error, actual in zip(errors, observed) if actual != 0]
    mean_observed = sum(observed) / len(observed)
    total_variation = sum((value - mean_observed) ** 2 for value in observed)
    r_squared = (
        1.0 - sum(error * error for error in errors) / total_variation
        if total_variation else 0.0
    )

    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=values.__getitem__)
        result = [0.0] * len(values)
        start = 0
        while start < len(order):
            end = start + 1
            while end < len(order) and values[order[end]] == values[order[start]]:
                end += 1
            rank = (start + end - 1) / 2
            for index in order[start:end]:
                result[index] = rank
            start = end
        return result

    left, right = ranks(observed), ranks(predicted)
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    denominator = math.sqrt(
        sum((a - left_mean) ** 2 for a in left)
        * sum((b - right_mean) ** 2 for b in right)
    )
    return {
        "count": float(len(rows)),
        "mae": mae,
        "rmse": rmse,
        "mape": sum(mape_values) / len(mape_values) if mape_values else 0.0,
        "r_squared": r_squared,
        "spearman": numerator / denominator if denominator else 0.0,
    }


def _fit_logistic(
    rows: list[dict[str, Any]],
    *,
    label_key: str = "execution_success",
    learning_rate: float = 0.03,
    iterations: int = 800,
) -> LogisticModel | None:
    labels = [_binary_label(row, label_key) for row in rows]
    if not rows or len(set(labels)) < 2:
        return None
    width = len(feature_vector(rows[0]))
    weights = [0.0] * width
    intercept = 0.0
    for _ in range(iterations):
        gradient = [0.0] * width
        intercept_gradient = 0.0
        for row, label in zip(rows, labels):
            values = feature_vector(row)
            raw = max(-30.0, min(30.0, intercept + sum(a * b for a, b in zip(weights, values))))
            probability = 1.0 / (1.0 + math.exp(-raw))
            error = probability - label
            intercept_gradient += error
            for index, value in enumerate(values):
                gradient[index] += error * value
        scale = learning_rate / len(rows)
        intercept -= scale * intercept_gradient
        weights = [weight - scale * value for weight, value in zip(weights, gradient)]
    return LogisticModel(tuple(weights), intercept, feature_names())


def _classification_metrics(
    model: LogisticModel | None,
    rows: list[dict[str, Any]],
    *,
    label_key: str = "execution_success",
) -> dict[str, float]:
    if model is None or not rows:
        return {"count": float(len(rows))}
    labels = [_binary_label(row, label_key) for row in rows]
    probabilities = [model.predict_probability(row) for row in rows]
    predicted = [probability >= 0.5 for probability in probabilities]
    true_positive = sum(choice and bool(label) for choice, label in zip(predicted, labels))
    true_negative = sum(not choice and not bool(label) for choice, label in zip(predicted, labels))
    false_positive = sum(choice and not bool(label) for choice, label in zip(predicted, labels))
    false_negative = sum(not choice and bool(label) for choice, label in zip(predicted, labels))
    positives = true_positive + false_negative
    negatives = true_negative + false_positive
    recall = true_positive / positives if positives else 0.0
    specificity = true_negative / negatives if negatives else 0.0
    return {
        "count": float(len(rows)),
        "positive_count": float(positives),
        "negative_count": float(negatives),
        "accuracy": (true_positive + true_negative) / len(rows),
        "precision": (
            true_positive / (true_positive + false_positive)
            if true_positive + false_positive else 0.0
        ),
        "recall": recall,
        "specificity": specificity,
        "balanced_accuracy": (
            (recall + specificity) / 2 if positives and negatives else 0.0
        ),
        "single_class_test": float(not positives or not negatives),
        "brier_score": sum((probability - label) ** 2
                           for probability, label in zip(probabilities, labels)) / len(rows),
        "log_loss": -sum(
            label * math.log(max(1e-15, min(1 - 1e-15, probability)))
            + (1 - label)
            * math.log(max(1e-15, min(1 - 1e-15, 1 - probability)))
            for probability, label in zip(probabilities, labels)
        ) / len(rows),
    }


def _binary_label(row: dict[str, Any], label_key: str) -> float:
    raw = row.get(label_key)
    if raw not in {None, ""}:
        if isinstance(raw, str):
            return float(raw.strip().lower() in {"1", "true", "yes", "ok"})
        return float(bool(raw))
    if label_key == "execution_success":
        return float(str(row.get("status", "")) == "ok")
    raise ValueError(f"row lacks binary label {label_key!r}")


def grouped_split(
    rows: list[dict[str, Any]],
    group_key: str,
    test_fraction: float,
    seed: int,
    held_out_groups: tuple[str, ...] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    groups = sorted({str(row.get(group_key, "")) for row in rows})
    if held_out_groups is not None:
        requested = set(held_out_groups)
        absent = requested - set(groups)
        if absent:
            raise ValueError(
                f"held-out groups are absent for {group_key}: {sorted(absent)}"
            )
        test_groups = requested
    else:
        random.Random(seed).shuffle(groups)
        test_count = max(1, round(len(groups) * test_fraction)) if len(groups) > 1 else 0
        test_groups = set(groups[:test_count])
    return (
        [row for row in rows if str(row.get(group_key, "")) not in test_groups],
        [row for row in rows if str(row.get(group_key, "")) in test_groups],
        sorted(test_groups),
    )


def train_cost_models(
    rows: list[dict[str, Any]],
    *,
    group_key: str = "workload_group",
    test_fraction: float = 0.2,
    seed: int = 20260728,
    alpha: float = 1.0,
    held_out_groups: tuple[str, ...] | None = None,
) -> CostModelBundle:
    rows = measured_training_rows(rows)
    if not rows:
        raise ValueError("cost-model training requires measured benchmark records")
    successful = [row for row in rows if str(row.get("status", "ok")) == "ok"]
    train, test, test_groups = grouped_split(
        successful, group_key, test_fraction, seed, held_out_groups
    )
    if not train:
        raise ValueError("cost-model training requires at least one training row")
    regressors: dict[str, LinearModel] = {}
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
        model = _fit_ridge(target_train, target, alpha)
        regressors[target] = model
        metrics[target] = _metrics(model, target_test or target_train, target)
    if not regressors:
        raise ValueError("dataset has no measurable cost-model targets")
    all_train, all_test, _ = grouped_split(
        rows, group_key, test_fraction, seed, held_out_groups
    )
    classifier = _fit_logistic(all_train)
    metrics["execution_success"] = _classification_metrics(classifier, all_test or all_train)
    error_rows = [
        row for row in rows
        if usable_binary_label(row, "error_constraint_feasible")
    ]
    try:
        error_train, error_test, error_test_groups = grouped_split(
            error_rows, group_key, test_fraction, seed, held_out_groups
        )
    except ValueError:
        # The regression/success split remains valid even if no labeled error
        # record exists in a requested transfer-test group.
        error_train, error_test, error_test_groups = error_rows, [], []
    error_classifier = _fit_logistic(
        error_train, label_key="error_constraint_feasible"
    )
    metrics["error_constraint_feasibility"] = _classification_metrics(
        error_classifier,
        error_test or error_train,
        label_key="error_constraint_feasible",
    )
    return CostModelBundle(
        feature_names=feature_names(),
        regressors=regressors,
        success_classifier=classifier,
        error_feasibility_classifier=error_classifier,
        metrics=metrics,
        split={
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
        target_scopes=target_scopes,
    )
