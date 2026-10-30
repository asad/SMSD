# MCP stdio container

This recipe installs a supplied SMSD 8.1.0 Linux CPython 3.12 wheel and its
`mcp` extra, then starts the actual `smsd-mcp` console entry point. It does not
build a native extension, start a web server or contact a model provider.

Prepare `wheelhouse/` beside the Dockerfile with the wheel for the image's
architecture and all compatible MCP dependency wheels. Acquire and review
these inputs separately. The recipe installs with `--no-index`; it cannot
fetch missing packages. Optional chemical standardisation additionally needs
RDKit and its dependencies, installed by a separately reviewed image recipe.
The default native representation tool needs no RDKit. The build runs
`copy_licences.py`, also beside the Dockerfile, which copies the `LICENSE`,
`NOTICE` and third-party licence files of the installed wheel to
`/usr/share/doc/smsd` in the image.

```sh
docker build --network=none --build-arg SMSD_WHEEL=YOUR_MATCHING_SMSD_WHEEL.whl \
  -t smsd-mcp:8.1.0 packaging/mcp
docker run --rm -i --network=none --read-only --tmpfs /tmp:rw,noexec,nosuid,size=64m \
  --cap-drop=ALL --security-opt=no-new-privileges --memory=1g --cpus=1 \
  smsd-mcp:8.1.0
```

Replace the wheel basename with the actual `smsd-8.1.0-...whl` file. Supply a
reviewed immutable `PYTHON_BASE` digest for a reproducible build. A cached base
image is required for an offline build; `--network=none` restricts build steps
and is not a guarantee that the Docker daemon cannot fetch a base image.
The runtime command uses standard input/output for an MCP client. Do not add
`-t`, which allocates a terminal instead of preserving the protocol stream.

The process runs as an unprivileged UID with no source or host-directory
mounts. The resource caps are example deployment limits; tool search budgets
apply after parsing and do not replace process limits. This source recipe is
not a published image or a claim that every base/architecture has been tested.
See [tool contracts and limits](../../docs/MCP_SERVER.md).
