#version 150

uniform mat4 modelViewProjectionMatrix;
in vec4 position;
//varying vec3 N; // normal

void main()
{
    //N = normalize(gl_Normal);
    gl_Position = modelViewProjectionMatrix * position;
}
