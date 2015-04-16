#version 150

uniform vec3 cameraPosition;

// Set by OpenFrameworks
uniform vec4 globalColor;
out vec4 outputColor;

void main()
{
    outputColor = globalColor;
}
