#version 150

uniform vec3 cameraPosition;

// Set by OpenFrameworks
uniform vec4 globalColor;
in vec4  color;
out vec4 outputColor;

/*
void main()
{
    //outputColor = globalColor;
    outputColor = color;
}
*/

//uniform vec3 Color;
//uniform vec3 lightDir;

void main(void)
{
    vec3 lightDir = vec3(10.0, -10.0, 10.0);
    lightDir = normalize(lightDir);
    
    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_PointCoord * 2.0 - vec2(1.0);
    
    float mag = dot(N.xy, N.xy);

    if (mag > 1.0) {
        discard; // kill pixels outside circle
    }
    
    N.z = sqrt(1.0 - mag);
    
    // calculate lighting
    float diffuse = max(0.0, dot(lightDir, N));
    
    outputColor = globalColor * diffuse;
}