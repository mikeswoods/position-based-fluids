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
#include "Definitions.h"

/******************************************************************************/

static float readParticleRadius(float currentSize);

/******************************************************************************/

using namespace std;

/*******************************************************************************
 * Simulation state
 ******************************************************************************/

void ofApp::setup()
{
    this->paused = true;
    this->advanceStep = false;
    
	ofSetLogLevel(OF_LOG_VERBOSE);
    ofSetVerticalSync(true);
    
    // This uses depth information for occlusion
    // rather than always drawing things on top of each other
    ofEnableDepthTest();
    
    // This sets the camera's distance from the object
    this->camera.setDistance(25);
    
    // Initialize from GL world:
    this->openCL.setupFromOpenGL();
    
    // Set the bounds of the simulation:
    AABB bounds(float3(-2.0f, -10.0f, -2.0f), float3(2.0f, 10.0f, 2.0f));

    // Instantiate the simulator:
    this->simulation =
        std::shared_ptr<Simulation>(new Simulation(this->openCL, bounds));
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
void ofApp::drawHeadsUpDisplay()
{
    // Show the current frame rate and frame count
    ofSetColor(255);
    ofFill();

    int hOffset     = 10;
    int vSpacing    = 15;
    int textYOffset = vSpacing;
    
    // FPS
    string fpsText = ofToString(ofGetFrameRate()) + " fps";
    ofDrawBitmapString(fpsText, hOffset, textYOffset += vSpacing);

    // Current frame
    ofDrawBitmapString("Frame: " + ofToString(this->simulation->getFrameNumber()), hOffset, textYOffset += vSpacing);

    // Hotkeys:
    ofDrawBitmapString("Hotkeys:", hOffset, textYOffset += vSpacing);
    vector<string> hotkeys;
    hotkeys.push_back("'s' = step");
    hotkeys.push_back("'z' = set particle radius");
    hotkeys.push_back("'p' or space = toggle pause");
    hotkeys.push_back("'r' = reset");
    hotkeys.push_back("'g' = toggle grid");
    hotkeys.push_back("'d' = toggle visual debugging");
    
    for (auto i = hotkeys.begin(); i != hotkeys.end(); i++) {
        ofDrawBitmapString(*i, hOffset, textYOffset += vSpacing);
    }
    
    // Paused flag
    if (this->isPaused()) {
        ofDrawBitmapString("Paused", hOffset, textYOffset += vSpacing);
    }
}

/**
 * Renders the simulation
 */
void ofApp::draw()
{
    ofBackground(0);
    
    this->camera.begin();
    this->simulation->draw();
    this->camera.end();
    
    this->drawHeadsUpDisplay();
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
        // Set particle radius
        case 'z':
            {
                float r = readParticleRadius(this->simulation->getParticleRadius());
                this->simulation->setParticleRadius(r);
                this->simulation->reset();
            }
            break;
        // Reset
        case 'r':
            {
                this->simulation->reset();
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
    string maxSizeStr     = ofToString(Constants::MAX_PARTICLE_MASS);
    string currentSizeStr = ofToString(currentSize);
    
    string input = ofSystemTextBoxDialog("Particle size? (> 0.0], max = " + maxSizeStr + ")", currentSizeStr);
    float r = ofToFloat(input);
    r = std::min(std::max(0.0f, r), Constants::MAX_PARTICLE_MASS);
    return r == 0.0f ? Constants::DEFAULT_PARTICLE_RADIUS : r;
}

/******************************************************************************/
