#version 150

uniform vec3 cameraPosition;

// Set by OpenFrameworks
uniform vec4 globalColor;
in vec4  color;
out vec4 outputColor;

void main()
{
    //outputColor = globalColor;
    outputColor = color;
}
