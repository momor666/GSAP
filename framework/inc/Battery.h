/**  Battery - Header
 *   @file       Battery.h
 *   @ingroup    GSAP-Support
 *
 *   @brief      Battery model class for prognostics
 *
 *   @author     Matthew Daigle
 *   @version    0.1.0
 *
 *   @pre        N/A
 *
 *      Contact: Matthew Daigle (matthew.j.daigle@nasa.gov)
 *      Created: March 5, 2016
 *
 *   @copyright Copyright (c) 2018 United States Government as represented by
 *     the Administrator of the National Aeronautics and Space Administration.
 *     All Rights Reserved.
 */

#ifndef BATTERY_H
#define BATTERY_H

#include <cmath>
#include <vector>

#include "ConfigMap.h"
#include "ModelFactory.h"
#include "PrognosticsModel.h"

// Default parameter values
static const double QMOBILE_DEFAULT_VALUE = 7600;

class Battery final : public PCOE::PrognosticsModel {
public:
    // Constructor
    Battery();

    // Constructor based on configMap
    Battery(const PCOE::ConfigMap& paramMap);

    // State indices
    struct stateIndices {
        static const unsigned int Tb  = 0;
        static const unsigned int Vo  = 1;
        static const unsigned int Vsn = 2;
        static const unsigned int Vsp = 3;
        static const unsigned int qnB = 4;
        static const unsigned int qnS = 5;
        static const unsigned int qpB = 6;
        static const unsigned int qpS = 7;
    };
    // Input indices
    struct inputIndices {
        static const unsigned int P = 0;
    };
    // Output indices
    struct outputIndices {
        static const unsigned int Tbm = 0;
        static const unsigned int Vm  = 1;
    };
    // Indices
    struct allIndices {
        struct stateIndices states;
        struct inputIndices inputs;
        struct outputIndices outputs;
    } indices;

    // Parameters
    struct Parameters {
        double An2;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double qnBMax;
        double U0p;
        double An7;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double Ro;
        //Internal Ohmic resistance of the battery.
        double Vol;
        double qnSMax;
        double F;
        double to;
        //Time constant for Ohmic potential.
        double Ap0;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double Ap9;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double An5;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double An9;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double qpBMax;
        double alpha;
        double VolSFraction;
        //Fraction of total electrode volume occupied by the surface control volume.
        double VEOD;
        //The voltage level that defines end-of-discharge (EOD).
        double qMax;
        double xpMin;
        //Minimum mole fraction for positive electrode.
        double Ap1;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double qpSMin;
        double An4;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double Ap3;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double qpSMax;
        double Ap4;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double An11;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double tsp;
        //Time constant for surface overpotential for positive electrode.
        double kn;
        //Lumped constant for Butler-Volmer equation for negative electrode.
        double Ap11;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double Ap5;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double kp;
        //Lumped constant for Butler-Volmer equation for positive electrode.
        double R;
        double qnBMin;
        double An12;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double An10;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double VolS;
        double xpMax;
        //Maximum mole fraction for positive electrode.
        double qBMax;
        double qSMax;
        double Ap8;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double An6;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double Sn;
        //Surface area for negative electrode.
        double qpMin;
        double Ap2;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double tsn;
        //Time constant for surface overpotential for negative electrode.
        double qnMin;
        double qpMax;
        double qnMax;
        double qnSMin;
        double U0n;
        //Empirical parameter for negative electrode in equilibrium potential.
        double qpBMin;
        double VolB;
        double Ap6;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double Sp;
        //Surface area for positive electrode.
        double Ap10;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double An0;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double xnMax;
        //Maximum mole fraction for negative electrode.
        double Ap12;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double An1;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double Ap7;
        //Empirical parameter for positive electrode in Redlich-Kister expansion.
        double An8;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
        double xnMin;
        //Minimum mole fraction for negative electrode.
        double tDiffusion;
        //Diffusion time constant (increasing this causes a decrease in diffusion rate).
        double qMobile;
        double An3;
        //Empirical parameter for negative electrode in Redlich-Kister expansion.
    } parameters;

    /** @brief      Execute state equation. This version of the function uses a given sampling time.
     *   @param      t Time
     *   @param      x Current state vector. This gets updated to the state at the new time.
     *   @param      u Input vector
     *   @param      n Process noise vector
     *   @param      dt Sampling time
     **/
    void stateEqn(const double t,
                  std::vector<double>& x,
                  const std::vector<double>& u,
                  const std::vector<double>& n,
                  const double dt);
    /** @brief      Execute output equation
     *   @param      t Time
     *   @param      x State vector
     *   @param      u Input vector
     *   @param      n Sensor noise vector
     *   @param      z Output vector. This gets updated to the new output at the given time.
     **/
    void outputEqn(const double t,
                   const std::vector<double>& x,
                   const std::vector<double>& u,
                   const std::vector<double>& n,
                   std::vector<double>& z);
    /** @brief      Execute threshold equation
     *   @param      t Time
     *   @param      x State vector
     *   @param      u Input vector
     **/
    bool thresholdEqn(const double t, const std::vector<double>& x, const std::vector<double>& u);
    /** @brief      Execute input equation.
     *               Determines what input (u) should be at the given time for the given input
     *parameters.
     *   @param      t Time
     *   @param      inputParameters Vector of input parameters, which are values that specify how
     *to define u for the given time.
     *   @param      u Input vector. Gets overwritten.
     **/
    void inputEqn(const double t,
                  const std::vector<double>& inputParameters,
                  std::vector<double>& u);
    /** @brief      Execute predicted output equation.
     *               Predicted outputs are those that are not measured, but are interested in being
     *predicted for prognostics.
     *   @param      t Time
     *   @param      x State vector
     *   @param      u Input vector
     *   @param      z Predicted output vector. Gets overwritten.
     **/
    void predictedOutputEqn(const double t,
                            const std::vector<double>& x,
                            const std::vector<double>& u,
                            std::vector<double>& z);

    // Set default parameters, based on 18650 cells
    void setParameters(const double qMobile = QMOBILE_DEFAULT_VALUE, const double Vol = 2e-5);

    /** @brief      Initialize state vector given initial inputs and outputs.
     *   @param      x Current state vector. This gets updated.
     *   @param      u Input vector
     *   @param      z Output vector
     **/
    void initialize(std::vector<double>& x,
                    const std::vector<double>& u,
                    const std::vector<double>& z);
};
#endif
