/*******************************************************************************
 * main.cpp
 * - Position based fluids entry point
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/


#include "ofMain.h"
#include "ofApp.h"

/******************************************************************************/

int main()
{
    ofSetCurrentRenderer(ofGLProgrammableRenderer::TYPE);
    ofSetupOpenGL(1024, 768, OF_WINDOW);
	ofRunApp(new ofApp());
}

/******************************************************************************/
