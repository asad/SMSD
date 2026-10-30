ARG PYTHON_BASE=python:3.12-slim
FROM ${PYTHON_BASE}

ARG SMSD_WHEEL
WORKDIR /app
COPY wheelhouse/ /opt/wheelhouse/
COPY requirements.txt app.py ./
RUN case "${SMSD_WHEEL}" in smsd-8.1.0-*.whl) ;; *) echo 'Supply a matching SMSD 8.1.0 CPython 3.12 Linux wheel' >&2; exit 1 ;; esac \
    && test "${SMSD_WHEEL##*/}" = "${SMSD_WHEEL}" \
    && python -m pip install --no-cache-dir --no-index --find-links=/opt/wheelhouse \
       "/opt/wheelhouse/${SMSD_WHEEL}" -r requirements.txt \
    && python -c "from importlib.metadata import version; assert version('smsd') == '8.1.0'" \
    && rm -rf /opt/wheelhouse

ENV HOME=/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 GRADIO_ANALYTICS_ENABLED=False \
    GRADIO_SERVER_NAME=0.0.0.0 GRADIO_SERVER_PORT=7860 OMP_NUM_THREADS=1
USER 10001:10001
EXPOSE 7860
CMD ["python", "app.py"]
