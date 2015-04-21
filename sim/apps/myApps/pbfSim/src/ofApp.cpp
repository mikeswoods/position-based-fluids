/*******************************************************************************
 * ofApp
 * - The OpenFrameworks application entry point
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "ofApp.h"
#include "Constants.h"

/******************************************************************************/

static float readParticleRadius(float currentSize);

/******************************************************************************/

using namespace std;

/*******************************************************************************
 * Simulation state
 ******************************************************************************/

void ofApp::setup()
{
    this->paused      = true;
    this->advanceStep = false;
    
#ifdef ENABLE_LOGGING
	ofSetLogLevel(OF_LOG_VERBOSE);
#endif
    
    // This uses depth information for occlusion
    // rather than always drawing things on top of each other
    ofEnableDepthTest();
    
    // Disable alpha blending, since we won't need it:
    ofDisableAlphaBlending();
    
    // This sets the camera's distance from the object
    if (!this->cameraSet) {
        this->camera.setDistance(50);
        this->cameraSet = true;
    }

    // Initialize from GL world:
    this->openCL.setupFromOpenGL();

#ifdef SIMPLE_SCENE

    AABB bounds(ofVec3f(-2.0f, 0.0f, -2.0f), ofVec3f(2.0f, 10.0f, 2.0f));

    int numParticles = 150;
    Parameters parameters = Constants::FOR_RADIUS_0_25;

    ofSetVerticalSync(true);

#else

    AABB bounds(ofVec3f(-30.0f, 0.0f, -10.0f), ofVec3f(30.0f, 50.0f, 10.0f));

    int numParticles = 10000;
    Parameters parameters = Constants::FOR_RADIUS_0_5;

    ofSetVerticalSync(false);

#endif

    ofLogNotice() << "sizeof(Parameters) = " << sizeof(Parameters) << endl
                  << parameters << endl;
    
    // Instantiate the simulator:
    this->simulation = std::shared_ptr<Simulation>(new Simulation(this->openCL
                                                                 ,bounds
                                                                 ,numParticles
                                                                 ,parameters));
}

/**
 * Called to update the state of the simulation
 */
void ofApp::update()
{
    if (this->isPaused()) {
        if (this->advanceStep) {
            this->simulation->step();
            this->advanceStep = false;
        }
    } else {
        this->simulation->step();
    }
}

/*******************************************************************************
 * Drawing
 ******************************************************************************/

/**
 * Draws a "heads up display" that shows the status of the simulation, as well
 * as some other pieces of pertinent information
 */
void ofApp::drawHeadsUpDisplay(ofEasyCam& camera) const
{
    ofVec3f cameraPos = camera.getPosition();
    ofVec3f targetPos = camera.getTarget().getPosition();
    
    // Show the current frame rate and frame count
    ofSetColor(255);
    ofFill();

    int hOffset     = 10;
    int vSpacing    = 15;
    int textYOffset = vSpacing;
    
    // Paused flag
    if (this->isPaused()) {
        ofSetColor(0, 255, 0);
        ofFill();
        ofDrawBitmapString("Paused", hOffset, textYOffset += vSpacing);
    }
    
    ofSetColor(255);
    ofFill();
    
    // FPS
    string fpsText = ofToString(ofGetFrameRate()) + " fps";
    ofDrawBitmapString(fpsText, hOffset, textYOffset += vSpacing);

    // Current frame
    ofDrawBitmapString("Frame: " + ofToString(this->simulation->getFrameNumber()), hOffset, textYOffset += vSpacing);

    // Hotkeys:
    ofDrawBitmapString("Hotkeys:", hOffset, textYOffset += vSpacing);
    vector<string> hotkeys;
    hotkeys.push_back("'s' = step");
    hotkeys.push_back("'p' or space = toggle pause");
    hotkeys.push_back("'r' = reset");
    hotkeys.push_back("'g' = toggle grid");
    hotkeys.push_back("'d' = toggle visual debugging");
    
    for (auto i = hotkeys.begin(); i != hotkeys.end(); i++) {
        ofDrawBitmapString(*i, hOffset, textYOffset += vSpacing);
    }

    // Status information:
    auto bounds = this->simulation->getBounds();
    auto minExt = bounds.getMinExtent();
    auto maxExt = bounds.getMaxExtent();
    
    // Camera Position:
    ofDrawBitmapString("Camera position: <" + ofToString(cameraPos.x) +
                                        "," + ofToString(cameraPos.y) +
                                        "," + ofToString(cameraPos.z) + ">"
                      ,hOffset
                      ,textYOffset += vSpacing);

    // Camera target:
    ofDrawBitmapString("Camera target: <" + ofToString(targetPos.x) +
                                      "," + ofToString(targetPos.y) +
                                      "," + ofToString(targetPos.z) + ">"
                       ,hOffset
                       ,textYOffset += vSpacing);

    // Bounds:
    ofDrawBitmapString("Bounds: <" +
                       ofToString(minExt.x) + "," + ofToString(minExt.y) + "," + ofToString(minExt.z) + "> <" +
                       ofToString(maxExt.x) + "," + ofToString(maxExt.y) + "," + ofToString(maxExt.z) + ">"
                      ,hOffset, textYOffset += vSpacing);
    // Cell count:
    auto cellsPerAxis = this->simulation->getCellsPerAxis();
    ofDrawBitmapString("Cells per axis: <" +
                       ofToString(cellsPerAxis[0]) + "," + ofToString(cellsPerAxis[1]) + "," + ofToString(cellsPerAxis[2]) + ">"
                       ,hOffset, textYOffset += vSpacing);
    
    // Particle count:
    ofDrawBitmapString("Particles: " + ofToString(this->simulation->getNumberOfParticles())
                      ,hOffset, textYOffset += vSpacing);
}

/**
 * Renders the simulation
 */
void ofApp::draw()
{
    ofBackground(0);
    
    this->camera.begin();

        this->simulation->draw(this->camera);

    this->camera.end();
    
    this->drawHeadsUpDisplay(this->camera);
}

/*******************************************************************************
 * Input handling
 ******************************************************************************/

/**
 * Tests if the simulation is paused
 */
bool ofApp::isPaused() const
{
    return this->paused;
}

/**
 * Enables/disabled pause state
 */
void ofApp::togglePaused()
{
    this->paused = !this->paused;
}

/**
 * Invokes a keypress callback function
 */
void ofApp::keyPressed(int key)
{
    switch (key) {
        // Animate bounds
        case 'a':
            {
                this->simulation->toggleBoundsAnimation();
            }
            break;
        // Pause
        case 'p':
        case ' ':
            {
                this->togglePaused();
            }
        break;
        // Step
        case 's':
            {
                this->advanceStep = true;
            }
        break;
        // Reset
        case 'r':
            {
                auto oldCam = this->camera;
                this->setup();
            }
            break;
        // Toggle visual debugging:
        case 'd':
            {
                this->simulation->toggleVisualDebugging();
            }
            break;
        // Toggle drawing the spatial grid:
        case 'g':
            {
                this->simulation->toggleDrawGrid();
            }
            break;
    }
}

void ofApp::keyReleased(int key)
{

}

void ofApp::mouseMoved(int x, int y )
{

}

void ofApp::mouseDragged(int x, int y, int button)
{

}

void ofApp::mousePressed(int x, int y, int button)
{

}

void ofApp::mouseReleased(int x, int y, int button)
{

}

void ofApp::windowResized(int w, int h)
{

}

void ofApp::gotMessage(ofMessage msg)
{

}

void ofApp::dragEvent(ofDragInfo dragInfo)
{

}

/*******************************************************************************
 * Utility functions
 ******************************************************************************/

float readParticleRadius(float currentSize)
{
    string maxSizeStr     = ofToString(Constants::MAX_PARTICLE_RADIUS);
    string currentSizeStr = ofToString(currentSize);
    
    string input = ofSystemTextBoxDialog("Particle size? (> 0.0], max = " + maxSizeStr + ")", currentSizeStr);
    float r = ofToFloat(input);
    r = std::min(std::max(0.0f, r), Constants::MAX_PARTICLE_RADIUS);
    return r == 0.0f ? Constants::DEFAULT_PARTICLE_RADIUS : r;
}

/******************************************************************************/
