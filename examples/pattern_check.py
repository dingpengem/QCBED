import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

pattern_file = Path("reference_data/pattern.dat")
output_file = Path("pattern_dat_inspection.png")

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------
pattern = np.loadtxt(pattern_file)

print("shape:", pattern.shape)
print("dtype:", pattern.dtype)
print("min:", np.min(pattern))
print("max:", np.max(pattern))
print("mean:", np.mean(pattern))
print("median:", np.median(pattern))

if pattern.shape != (256, 256):
    print("Warning: pattern shape is not 256 x 256")

ny, nx = pattern.shape

# ------------------------------------------------------------
# Display versions
# ------------------------------------------------------------

# Linear
linear = pattern

# Log display only
log_display = np.log1p(pattern - np.min(pattern))

# Percentile contrast
vmin, vmax = np.percentile(pattern, [1, 99.5])

# ------------------------------------------------------------
# Central profiles
# ------------------------------------------------------------
cx = nx // 2
cy = ny // 2

profile_x = pattern[cy, :]
profile_y = pattern[:, cx]

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Linear
im0 = axes[0, 0].imshow(
    linear,
    cmap="gray",
    origin="upper",
)
axes[0, 0].axvline(cx)
axes[0, 0].axhline(cy)
axes[0, 0].set_title("Linear intensity")
axes[0, 0].set_xlabel("x (pixel)")
axes[0, 0].set_ylabel("y (pixel)")
fig.colorbar(im0, ax=axes[0, 0])

# Percentile
im1 = axes[0, 1].imshow(
    pattern,
    cmap="gray",
    origin="upper",
    vmin=vmin,
    vmax=vmax,
)
axes[0, 1].axvline(cx)
axes[0, 1].axhline(cy)
axes[0, 1].set_title("1–99.5 percentile contrast")
axes[0, 1].set_xlabel("x (pixel)")
axes[0, 1].set_ylabel("y (pixel)")
fig.colorbar(im1, ax=axes[0, 1])

# Log
im2 = axes[1, 0].imshow(
    log_display,
    cmap="gray",
    origin="upper",
)
axes[1, 0].axvline(cx)
axes[1, 0].axhline(cy)
axes[1, 0].set_title("Log display")
axes[1, 0].set_xlabel("x (pixel)")
axes[1, 0].set_ylabel("y (pixel)")
fig.colorbar(im2, ax=axes[1, 0])

# Center profiles
axes[1, 1].plot(np.arange(nx), profile_x, label=f"horizontal, y={cy}")
axes[1, 1].plot(np.arange(ny), profile_y, label=f"vertical, x={cx}")
axes[1, 1].set_title("Central intensity profiles")
axes[1, 1].set_xlabel("pixel")
axes[1, 1].set_ylabel("intensity")
axes[1, 1].legend()

plt.tight_layout()
plt.savefig(output_file, dpi=200)
plt.show()

print(f"Saved: {output_file}")