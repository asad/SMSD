ARG PYTHON_BASE=python:3.12-slim
FROM ${PYTHON_BASE}

ARG SMSD_WHEEL
COPY wheelhouse/ /opt/wheelhouse/
COPY copy_licences.py /opt/copy_licences.py
RUN case "${SMSD_WHEEL}" in smsd-8.1.0-*.whl) ;; *) echo 'Set SMSD_WHEEL to the matching SMSD 8.1.0 wheel basename' >&2; exit 1 ;; esac \
    && test "${SMSD_WHEEL##*/}" = "${SMSD_WHEEL}" \
    && python -m pip install --no-cache-dir --no-index --find-links=/opt/wheelhouse \
       "/opt/wheelhouse/${SMSD_WHEEL}[mcp]" \
    && python -c "from importlib.metadata import version; assert version('smsd') == '8.1.0'; assert version('mcp')" \
    && python /opt/copy_licences.py /usr/share/doc/smsd \
    && rm -rf /opt/wheelhouse /opt/copy_licences.py

ENV HOME=/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1
USER 10001:10001
ENTRYPOINT ["smsd-mcp"]
