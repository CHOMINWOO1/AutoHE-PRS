FROM python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir openfhe==1.5.1.0.22.4

WORKDIR /opt/autohe-prs
COPY pyproject.toml README.md ./
COPY src ./src
RUN python -m pip install --no-cache-dir .

WORKDIR /work
ENTRYPOINT ["autohe-prs"]
CMD ["--help"]
