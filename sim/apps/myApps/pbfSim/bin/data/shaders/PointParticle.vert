#version 150

uniform float particleRadius;
uniform vec3  cameraPosition;
uniform mat4  modelViewProjectionMatrix;

out vec4 color;
in vec4 position;
in vec4 normal;

// Constants (tweakable):
const float minPointScale = 0.01;
const float maxPointScale = 25.0;
const float maxDistance   = 100.0;

void main()
{
    gl_Position = modelViewProjectionMatrix * position;

    float cameraDist = distance(gl_Position.xyz, cameraPosition);
    float pointScale = 1.0 - clamp(cameraDist / maxDistance, 0.001, 1.0);
    pointScale = max(pointScale, minPointScale);
    pointScale = min(pointScale, maxPointScale);
    
    float wScale = position.w;
    if (wScale == 0.0) {
        wScale = 0.0001;
    }
    
    gl_PointSize = (wScale * pointScale) * particleRadius;
    
    color = vec4(1.0, 1.0, 1.0, 1.0);
}
