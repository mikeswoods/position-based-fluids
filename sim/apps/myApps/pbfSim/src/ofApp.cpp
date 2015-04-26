/*******************************************************************************
 * ofApp.cpp
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

using namespace std;

/*******************************************************************************
 * Simulation state
 ******************************************************************************/

void ofApp::setup()
{
#ifdef ENABLE_LOGGING
    ofSetLogLevel(OF_LOG_VERBOSE);
#endif

    ofSetVerticalSync(true);
    
    // disable alpha blending, since we won't need it:

    ofDisableAlphaBlending();
    
    //ofEnableDepthTest();
    
    // set up the UI controls:
    
    this->gui.setup("Animation");
    this->gui.add(this->toggleAnimateBounds.setup("Toggle", false));
    this->gui.add(this->toggleAnimateBothSides.setup("2-sides", false));
    this->gui.add(this->selectSineAnim.setup("~ Sine wave"));
    this->gui.add(this->selectRampAnim.setup("~ Linear ramp"));
    this->gui.add(this->selectCompressAnim.setup("~ Compress"));
    this->gui.add(this->periodAnimSlider.setup("Period", 1.0f, 0.01f, 2.0f));
    this->gui.add(this->ampAnimSlider.setup("Amplitude", 1.0f, 1.0f, 15.0f));
    this->gui.add(this->resetBounds.setup("Reset Bounds"));

    this->selectSineAnim.addListener(this, &ofApp::setSineAnim);
    this->selectRampAnim.addListener(this, &ofApp::setRampAnim);
    this->selectCompressAnim.addListener(this, &ofApp::setCompressAnim);
    this->periodAnimSlider.addListener(this, &ofApp::setAnimPeriod);
    this->ampAnimSlider.addListener(this, &ofApp::setAnimAmp);
    this->resetBounds.addListener(this, &ofApp::doResetBounds);
    
    // and the flag initial values:
    
    this->paused      = true;
    this->advanceStep = false;
    
    // set the camera's distance from the object:

    if (!this->cameraSet) {
        this->camera.setDistance(50);
        this->cameraSet = true;
    }

    // Initialize from GL world:

    this->openCL.setupFromOpenGL();

    // Set the scene parameters:
    
#ifdef SIMPLE_SCENE

    AABB bounds(ofVec3f(-2.0f, 0.0f, -2.0f), ofVec3f(2.0f, 10.0f, 2.0f));
    int numParticles      = 47;
    Parameters parameters = Constants::DEFAULT_PARAMS;

    this->simulation = shared_ptr<Simulation>(new Simulation(this->openCL
                                                            ,bounds
                                                            ,numParticles
                                                            ,Constants::DEFAULT_DT
                                                            ,ofVec3f(2,2,2)
                                                            ,parameters));
    
#else

    AABB bounds(ofVec3f(-30.0f, -10.0f, -10.0f), ofVec3f(30.0f, 80.0f, 10.0f));
    int numParticles      = 10000;
    Parameters parameters = Constants::DEFAULT_PARAMS;

    this->simulation = shared_ptr<Simulation>(new Simulation(this->openCL
                                                            ,bounds
                                                            ,numParticles
                                                            ,parameters));
#endif
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

    // Animate bounds?

    if (this->toggleAnimateBounds) {
        this->simulation->enableBoundsAnimation();
    } else {
        this->simulation->disableBoundsAnimation();
    }
    
    if (this->toggleAnimateBothSides) {
        this->simulation->enableBothSidesAnimation();
    } else {
        this->simulation->disableBothSidesAnimation();
    }
}

/*******************************************************************************
 * Control callback functions
 ******************************************************************************/

void ofApp::setSineAnim()
{
    this->simulation->setAnimationType(Simulation::SINE_WAVE);
}

void ofApp::setRampAnim()
{
    this->simulation->setAnimationType(Simulation::LINEAR_RAMP);
}

void ofApp::setCompressAnim()
{
    this->simulation->setAnimationType(Simulation::COMPRESS);
}

void ofApp::doResetBounds()
{
    this->simulation->disableBoundsAnimation();
    this->simulation->resetBounds();
}

void ofApp::setAnimPeriod(float& period)
{
    this->simulation->setAnimationPeriod(period);
}

void ofApp::setAnimAmp(float& amp)
{
    this->simulation->setAnimationAmp(amp);
}

/*******************************************************************************
 * Drawing
 ******************************************************************************/

/**
 * Draws a "heads up display" that shows the status of the simulation, as well
 * as some other pieces of pertinent information
 */
void ofApp::drawHeadsUpDisplay(ofEasyCam& camera)
{
    // Show the controls:
    
    ofVec3f cameraPos = camera.getPosition();
    
    // Show the current frame rate and frame count

    int w           = static_cast<int>(ofGetScreenWidth());
    int hOffset     = static_cast<int>(static_cast<float>(w) - (0.425f * static_cast<float>(w)));
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

    // Hotkeys

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

    // Status information

    auto bounds = this->simulation->getBounds();
    auto minExt = bounds.getMinExtent();
    auto maxExt = bounds.getMaxExtent();
    
    // Camera Position

    ofDrawBitmapString("Camera position: <" + ofToString(cameraPos.x) +
                                        "," + ofToString(cameraPos.y) +
                                        "," + ofToString(cameraPos.z) + ">"
                      ,hOffset
                      ,textYOffset += vSpacing);

    // Bounds

    ofDrawBitmapString("Bounds: <" +
                       ofToString(minExt.x) + "," + ofToString(minExt.y) + "," + ofToString(minExt.z) + "> <" +
                       ofToString(maxExt.x) + "," + ofToString(maxExt.y) + "," + ofToString(maxExt.z) + ">"
                      ,hOffset, textYOffset += vSpacing);

    // Cell count

    auto cellsPerAxis = this->simulation->getCellsPerAxis();
    ofDrawBitmapString("Cells per axis: <" +
                       ofToString(cellsPerAxis[0]) + "," + ofToString(cellsPerAxis[1]) + "," + ofToString(cellsPerAxis[2]) + ">"
                       ,hOffset, textYOffset += vSpacing);
    
    // Particle count

    ofDrawBitmapString("Particles: " + ofToString(this->simulation->getNumberOfParticles())
                      ,hOffset, textYOffset += vSpacing);

    ofDisableDepthTest();
    this->gui.draw();
    ofEnableDepthTest();
}

/**
 * Renders the simulation
 */
void ofApp::draw()
{
    ofBackground(0);
    
    // Render the current step of the simulation:
    
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

/******************************************************************************/
