from colorsys import hls_to_rgb

def rainbow_color_steps(steps:int, end=6/8):
    raw_colors = [hls_to_rgb(end * i/(n-1), 0.5, 1) for i in range(steps)]
    return [[int(175 * x) for x in color] for color in raw_colors]


print (rainbow_color_stops(100))