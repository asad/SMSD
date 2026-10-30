# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Start the MCP server over standard input and output, list its tools and run one search."""

import asyncio
import json
import os
import sys

try:
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client
except ImportError:
    sys.exit('This example needs the MCP SDK: python -m pip install "smsd[mcp]"')


async def main():
    # The server runs with this interpreter and environment, so it imports the same smsd.
    server = StdioServerParameters(command=sys.executable, args=["-m", "smsd.mcp_server"],
                                   env=dict(os.environ))
    async with stdio_client(server) as (reader, writer):
        async with ClientSession(reader, writer) as session:
            await session.initialize()
            tools = await session.list_tools()
            print("Tools:", ", ".join(tool.name for tool in tools.tools))
            response = await session.call_tool("find_mcs", {"query_smiles": "CCO",
                                                            "target_smiles": "CCN"})
            result = json.loads(response.content[0].text)
            print(f"find_mcs CCO/CCN: status={result['status']}, optimal={result['optimal']}, "
                  f"complete={result['complete']}")


if __name__ == "__main__":
    asyncio.run(main())
