"""tjophot -- a photometry pipeline for TJO/MEIA2 supernova imaging.

The pipeline is deliberately instrument-agnostic at the edges: everything that
knows about TJO lives in ``config/instruments/tjo_meia2.yaml``.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
