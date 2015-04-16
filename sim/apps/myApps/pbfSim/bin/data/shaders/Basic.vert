#version 150

uniform mat4 modelViewProjectionMatrix;
out vec4 color;
in vec4 position;
in vec4 normal;
//in vec4 texcoord

void main()
{
    gl_Position = modelViewProjectionMatrix * position;
    color = normalize(normal);
}
