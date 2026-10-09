# Save a matplotlib figure as a journal-ready TIFF: 600 dpi, RGB (no alpha channel), LZW-compressed.
import io

from PIL import Image


def save_tiff(fig, path, dpi=600):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    Image.open(buf).convert("RGB").save(path, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
