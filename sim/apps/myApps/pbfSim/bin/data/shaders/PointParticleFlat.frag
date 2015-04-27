#version 150

uniform vec4 globalColor;
in vec4  color;
out vec4 outputColor;

const float OUTER_BORDER = 0.1;
const float INNER_BORDER = 0.4;

/**
 * Adapted from http://mmmovania.blogspot.de/2011/01/point-sprites-as-spheres-in-opengl33.html
 * Copyright (C) 2011 - Movania Muhammad Mobeen
 */
void main(void)
{
    vec3 N;
    N.xy = gl_PointCoord * 2.0 - vec2(1.0);
    
    float mag = 1.0 - dot(N.xy, N.xy);
    
    if (mag <= OUTER_BORDER) {
        discard;
    } else if (mag <= INNER_BORDER) {
        outputColor = vec4(0.58984375, 0.51171875, 0.62890625, 1.0);
    } else {
        outputColor = vec4(0.46484375, 0.80859375, 0.78515625, 1.0);
    }
}
