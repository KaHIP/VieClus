#!/usr/bin/env python3
"""Generate the GitHub social preview image (1280x640) for VieClus."""

from PIL import Image, ImageDraw, ImageFont
import math
import os
import random

random.seed(42)

W, H = 1280, 640
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "social-preview.png")

canvas = Image.new("RGB", (W, H), (26, 29, 35))
draw = ImageDraw.Draw(canvas)

# Subtle radial gradient background
cx, cy = W // 2, H // 2
for r in range(500, 0, -2):
    frac = r / 500
    c1 = int(26 + (40 - 26) * (1 - frac))
    c2 = int(29 + (45 - 29) * (1 - frac))
    c3 = int(35 + (55 - 35) * (1 - frac))
    draw.ellipse([cx - r * 1.6, cy - r, cx + r * 1.6, cy + r], fill=(c1, c2, c3))

# Cluster colors (matching VieClus logo: red, blue, yellow)
cluster_red = (231, 76, 60)
cluster_blue = (52, 152, 219)
cluster_yellow = (241, 196, 15)
inter_edge_color = (136, 136, 136)

# Graph vertices arranged in 3 clusters (left side of image)
# Red cluster (top-left)
red_nodes = [
    (100, 140), (180, 100), (260, 120), (200, 200), (120, 230),
    (280, 190), (160, 160),
]
# Blue cluster (top-right area of graph section)
blue_nodes = [
    (370, 90), (440, 130), (350, 170), (420, 210), (480, 170),
]
# Yellow cluster (bottom area)
yellow_nodes = [
    (200, 340), (280, 310), (350, 350), (300, 410), (220, 420),
    (380, 290), (420, 370),
]

all_nodes = []
node_cluster = []
for n in red_nodes:
    all_nodes.append(n)
    node_cluster.append(0)
for n in blue_nodes:
    all_nodes.append(n)
    node_cluster.append(1)
for n in yellow_nodes:
    all_nodes.append(n)
    node_cluster.append(2)

cluster_colors = [cluster_red, cluster_blue, cluster_yellow]

# Generate intra-cluster edges (nearby nodes within same cluster)
intra_edges = []
inter_edges = []
for i in range(len(all_nodes)):
    for j in range(i + 1, len(all_nodes)):
        x1, y1 = all_nodes[i]
        x2, y2 = all_nodes[j]
        dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        if node_cluster[i] == node_cluster[j]:
            if dist < 160:
                intra_edges.append((i, j))
        else:
            if dist < 120 and random.random() < 0.3:
                inter_edges.append((i, j))

# Draw inter-cluster edges (dashed-like, dim)
for i, j in inter_edges:
    x1, y1 = all_nodes[i]
    x2, y2 = all_nodes[j]
    # Draw dashed line
    length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    dash_len = 6
    gap_len = 6
    steps = int(length / (dash_len + gap_len))
    for s in range(steps):
        t1 = s * (dash_len + gap_len) / length
        t2 = min((s * (dash_len + gap_len) + dash_len) / length, 1.0)
        sx = x1 + (x2 - x1) * t1
        sy = y1 + (y2 - y1) * t1
        ex = x1 + (x2 - x1) * t2
        ey = y1 + (y2 - y1) * t2
        draw.line([(sx, sy), (ex, ey)], fill=(100, 100, 110, 80), width=1)

# Draw intra-cluster edges
for i, j in intra_edges:
    x1, y1 = all_nodes[i]
    x2, y2 = all_nodes[j]
    c = cluster_colors[node_cluster[i]]
    edge_c = tuple(int(v * 0.5) for v in c)
    draw.line([(x1, y1), (x2, y2)], fill=edge_c, width=2)

# Glow for all cluster nodes
for idx, (x, y) in enumerate(all_nodes):
    c = cluster_colors[node_cluster[idx]]
    bg = (26, 29, 35)
    for r in range(22, 5, -1):
        alpha_frac = 1 - (r - 5) / 17
        gc = tuple(int(c[k] * 0.15 * alpha_frac + (1 - 0.15 * alpha_frac) * bg[k]) for k in range(3))
        draw.ellipse([x - r, y - r, x + r, y + r], fill=gc)

# Draw nodes
for idx, (x, y) in enumerate(all_nodes):
    c = cluster_colors[node_cluster[idx]]
    bright = tuple(min(255, int(v * 1.3)) for v in c)
    r = 9
    draw.ellipse([x - r, y - r, x + r, y + r], fill=c, outline=bright, width=2)

# Fonts
try:
    font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 80)
    font_sub = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 28)
    font_tag = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
    font_legend = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 18)
except OSError:
    font_title = font_sub = font_tag = font_legend = ImageFont.load_default()

text_x = 600

# Title
draw.text((text_x, 160), "VieClus", fill=(240, 245, 255), font=font_title)

# Separator
draw.line([(text_x, 265), (text_x + 540, 265)], fill=(60, 80, 110), width=2)

# Subtitle
draw.text((text_x, 285), "Vienna Graph Clustering", fill=(180, 195, 220), font=font_sub)

# Tagline
draw.text((text_x, 350), "Multilevel memetic algorithm for", fill=(110, 130, 160), font=font_tag)
draw.text((text_x, 380), "graph clustering / community detection", fill=(110, 130, 160), font=font_tag)

# Legend
legend_y = 460
legend_items = [
    (cluster_red, "Cluster 1"),
    (cluster_blue, "Cluster 2"),
    (cluster_yellow, "Cluster 3"),
]
lx = text_x
for color, label in legend_items:
    draw.ellipse([lx - 5, legend_y - 5, lx + 5, legend_y + 5], fill=color)
    draw.text((lx + 14, legend_y - 10), label, fill=(136, 144, 160), font=font_legend)
    lx += 140

# Dashed line legend
dash_y = legend_y + 35
draw.line([(text_x, dash_y), (text_x + 30, dash_y)], fill=(136, 136, 136), width=1)
draw.text((text_x + 40, dash_y - 10), "Inter-cluster edges", fill=(136, 144, 160), font=font_legend)

# Decorative dots at bottom
for i, dx in enumerate(range(0, 340, 40)):
    dot_x = text_x + dx
    dot_y = 560
    colors = [cluster_red, cluster_blue, cluster_yellow]
    if i % 3 == 0:
        draw.ellipse([dot_x - 5, dot_y - 5, dot_x + 5, dot_y + 5], fill=colors[0])
    elif i % 3 == 1:
        draw.ellipse([dot_x - 4, dot_y - 4, dot_x + 4, dot_y + 4], fill=colors[1])
    else:
        draw.ellipse([dot_x - 4, dot_y - 4, dot_x + 4, dot_y + 4], fill=colors[2])

canvas.save(OUT_PATH, "PNG", quality=95)
print(f"Saved {OUT_PATH}")
