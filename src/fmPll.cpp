#include <cmath>
#include <vector>
#include "dy4.h"

std::vector<double> fmPll(const std::vector<double>& pllIn, double freq, double Fs, double ncoScale = 2.0, double phaseAdjust = 0.0, double normBandwidth = 0.01) {

    // scale factors for proportional/integrator terms
    // these scale factors were derived assuming the following:
    // damping factor of 0.707 (1 over square root of 2)
    // there is no oscillator gain and no phase detector gain
    double Cp = 2.666;
    double Ci = 3.555;

    // gain for the proportional term
    double Kp = normBandwidth * Cp;
    // gain for the integrator term
    double Ki = normBandwidth * normBandwidth * Ci;

    // output vector for the NCO
    std::vector<double> ncoOut(pllIn.size() + 1);

    // initialize internal state
    double integrator = 0.0;
    double phaseEst = 0.0;
    double feedbackI = 1.0;
    double feedbackQ = 0.0;
    ncoOut[0] = 1.0;
    int trigOffset = 0;
    // note: state saving will be needed for block processing

    for (int k = 0; k < pllIn.size(); ++k) {

        // phase detector
        double errorI = pllIn[k] * (+feedbackI); // complex conjugate of the
        double errorQ = pllIn[k] * (-feedbackQ); // feedback complex exponential

        // four-quadrant arctangent discriminator for phase error detection
        double errorD = std::atan2(errorQ, errorI);

        // loop filter
        integrator = integrator + Ki * errorD;

        // update phase estimate
        phaseEst = phaseEst + Kp * errorD + integrator;

        // internal oscillator
        ++trigOffset;
        double trigArg = 2 * PI * (freq / Fs) * (trigOffset) + phaseEst;
        feedbackI = std::cos(trigArg);
        feedbackQ = std::sin(trigArg);
        ncoOut[k + 1] = std::cos(trigArg * ncoScale + phaseAdjust);

    }

    // for stereo only the in-phase NCO component should be returned
    // for block processing you should also return the state
    return ncoOut;
    // for RDS add also the quadrature NCO component to the output

}

int main() {

    return 0;

}