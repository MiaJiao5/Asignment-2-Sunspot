# /// script
# requires-python = ">=3.10"
# dependencies = ["requests"]
# ///

"""
Fetch the numbers once, save the raw reply to data/, and never fetch again.

    uv run fetch.py

Change URL and FILE. The default is the Hong Kong Observatory's daily mean
temperature for 2026, so the template runs before you have touched it and you
can see what a file looks like when it arrives. It is an example, not your
phenomenon: handing it in unchanged is handing in nothing.
"""

from pathlib import Path

import requests

URL = "https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv"            # SILSO daily total sunspot numbers
FILE = "SN_d_tot_V2.0.csv"                                           # keep the publisher's filename
                                                                      # keep the publisher's extension
HERE = Path(__file__).parent
DATA = HERE / "data"


def fetch(url, path):
    """Ask for the file once. If it is already in data/, do nothing."""
    if path.exists():
        print(f"data/{path.name} is already here ({path.stat().st_size // 1024} KB). "
              "Delete it to fetch again.")
        return path
    DATA.mkdir(exist_ok=True)
    print(f"asking {url}")
    reply = requests.get(url, timeout=60, headers={"User-Agent": "SD5913 PolyU student"})
    reply.raise_for_status()
    path.write_bytes(reply.content)      # the raw reply, byte for byte: what arrived is what gets committed
    print(f"saved data/{path.name} ({path.stat().st_size // 1024} KB). Now: git add data")
    return path


if __name__ == "__main__":
    fetch(URL, DATA / FILE)
