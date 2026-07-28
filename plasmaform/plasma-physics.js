(function (global) {
  'use strict';

  // ── CODATA 2022 physical constants (SI) ──
  const constants = {
    elementaryCharge:   1.602176634e-19,   // C (exact)
    electronMass:       9.1093837139e-31,  // kg
    protonMass:         1.67262192595e-27, // kg
    vacuumPermeability: 1.25663706127e-6,  // N/A^2 (measured, CODATA 2022)
    vacuumPermittivity: 8.8541878188e-12,  // F/m (derived)
    speedOfLight:        2.99792458e8,      // m/s (exact)
    boltzmannConstant:  1.380649e-23,       // J/K (exact)
  };

  const { elementaryCharge: E, electronMass: ME, protonMass: MP,
          vacuumPermeability: MU0, vacuumPermittivity: EPS0,
          speedOfLight: C } = constants;

  // ── Unit conversions ──
  const conversions = {
    gaussToTesla:  G  => G * 1e-4,
    teslaToGauss:  T  => T * 1e4,
    perCcToPerM3:  n  => n * 1e6,
    perM3ToPerCc:  n  => n * 1e-6,
    eVToJoule:     eV => eV * E,
  };

  function requirePositiveFinite(value, label) {
    if (!Number.isFinite(value) || value <= 0) {
      throw new Error(`${label} must be a positive number.`);
    }
  }

  // ── Angular frequencies ──
  function electronGyroAngular(B) {
    // omega_ce = eB/m_e
    return (E * B) / ME;
  }

  function ionGyroAngular(B, Z, mu) {
    // omega_ci = ZeB/(mu m_p)
    return (Z * E * B) / (mu * MP);
  }

  function electronPlasmaAngular(ne) {
    // omega_pe = sqrt(ne e^2 / (eps0 m_e))
    return Math.sqrt((ne * E * E) / (EPS0 * ME));
  }

  function ionPlasmaAngular(ni, Z, mu) {
    // omega_pi = sqrt(ni Z^2 e^2 / (eps0 mu m_p))
    return Math.sqrt((ni * Z * Z * E * E) / (EPS0 * mu * MP));
  }

  function lowerHybridAngular(B, ne, Z, mu) {
    // omega_LH = sqrt[ omega_ci omega_ce / (1 + omega_ce^2/omega_pe^2) ]
    const wce = electronGyroAngular(B);
    const wci = ionGyroAngular(B, Z, mu);
    const wpe = electronPlasmaAngular(ne);
    return Math.sqrt((wci * wce) / (1 + (wce * wce) / (wpe * wpe)));
  }

  // ── Length scales ──
  function electronDebyeLength(TeEv, ne) {
    // lambda_De = sqrt( eps0 * Te[J] / (ne e) ),  Te[J] = TeEv * e
    return Math.sqrt((EPS0 * TeEv * E) / (ne * E * E));
  }

  function ionGyroradius(TiEv, B, Z, mu) {
    // rho_i = v_Ti / Omega_ci ,  v_Ti = sqrt(Ti[J]/m_i)  (NRL one-factor convention)
    const vTi = Math.sqrt((TiEv * E) / (mu * MP));
    const wci = ionGyroAngular(B, Z, mu);
    return vTi / wci;
  }

  function electronGyroradius(TeEv, B) {
    // r_e = v_Te / Omega_ce
    const vTe = electronThermalSpeed(TeEv);
    const wce = electronGyroAngular(B);
    return vTe / wce;
  }

  function ionInertialLength(ni, Z, mu) {
    // d_i = c / omega_pi
    return C / ionPlasmaAngular(ni, Z, mu);
  }

  function electronInertialLength(ne) {
    // d_e = c / omega_pe
    return C / electronPlasmaAngular(ne);
  }

  // ── Velocities ──
  function alfvenSpeed(B, n, mu) {
    // v_A = B / sqrt(mu0 n mu m_p)   (n = ion number density, mu = ion mass number)
    return B / Math.sqrt(MU0 * n * mu * MP);
  }

  function ionSoundSpeed(TeEv, Z, mu, gamma) {
    // C_s^IA = sqrt( gamma Z Te[J] / (mu m_p) )
    return Math.sqrt((gamma * Z * TeEv * E) / (mu * MP));
  }

  function mhdSoundSpeed(TeEv, TiEv, Z, mu, gammaElectron, gammaIon) {
    // C_s^MHD = sqrt( (gammaElectron*Z*Te[J] + gammaIon*Ti[J]) / (mu m_p) )
    return Math.sqrt(((gammaElectron * Z * TeEv + gammaIon * TiEv) * E) / (mu * MP));
  }

  function electronThermalSpeed(TeEv) {
    // v_Te = sqrt( Te[J] / m_e )   (NRL one-factor convention, not sqrt(2kT/m))
    return Math.sqrt((TeEv * E) / ME);
  }

  function ionThermalSpeed(TiEv, mu) {
    return Math.sqrt((TiEv * E) / (mu * MP));
  }

  // ── Dimensionless parameters ──
  function speciesBeta(n, TEv, B) {
    // beta_s = 2 mu0 n_s Te[J] / B^2
    return (2 * MU0 * n * TEv * E) / (B * B);
  }

  function electronIonMassRatio(mu) {
    return ME / (mu * MP);
  }

  function kawRegime(betaE, mu) {
    const ratio = betaE / electronIonMassRatio(mu);
    let label;
    if (!Number.isFinite(ratio)) label = 'Undefined';
    else if (ratio > 10) label = 'Kinetic Alfvén limit';
    else if (ratio < 0.1) label = 'Inertial Alfvén limit';
    else label = 'Transitional regime';
    return { ratio, label };
  }

  // ── KAW / Walén diagnostics ──
  function alfvenEquivalentVelocity(dB, n, mu) {
    // delta v_A = delta B / sqrt(mu0 rho),  rho = n mu m_p  (sign of dB preserved)
    return dB / Math.sqrt(MU0 * n * mu * MP);
  }

  // ── Master state calculator used by the Calculator / Validation / Worked-Example tabs ──
  function plasmaState(params) {
    const {
      niCm3, BnT, TeEv, TiEv, VswKms,
      Z = 1, mu = 1,
      gammaSound = 1, gammaElectron = 5 / 3, gammaIon = 5 / 3,
    } = params || {};

    requirePositiveFinite(niCm3, 'Ion density');
    requirePositiveFinite(BnT, 'Magnetic field');
    requirePositiveFinite(TeEv, 'Electron temperature');
    requirePositiveFinite(TiEv, 'Ion temperature');
    requirePositiveFinite(Z, 'Charge state');
    requirePositiveFinite(mu, 'Ion mass number');
    if (!Number.isFinite(VswKms)) {
      throw new Error('Solar wind speed must be a number.');
    }

    const B  = BnT * 1e-9; // nT -> T
    const ni = conversions.perCcToPerM3(niCm3);
    const ne = Z * ni;

    const wce = electronGyroAngular(B);
    const wci = ionGyroAngular(B, Z, mu);
    const wpe = electronPlasmaAngular(ne);
    const wpi = ionPlasmaAngular(ni, Z, mu);
    const wLH = lowerHybridAngular(B, ne, Z, mu);

    const fce = wce / (2 * Math.PI);
    const fci = wci / (2 * Math.PI);
    const fpe = wpe / (2 * Math.PI);
    const fpi = wpi / (2 * Math.PI);
    const fLH = wLH / (2 * Math.PI);

    const lDe = electronDebyeLength(TeEv, ne);
    const ri  = ionGyroradius(TiEv, B, Z, mu);
    const re  = electronGyroradius(TeEv, B);
    const di  = ionInertialLength(ni, Z, mu);
    const de  = electronInertialLength(ne);

    const vA    = alfvenSpeed(B, ni, mu);
    const csIA  = ionSoundSpeed(TeEv, Z, mu, gammaSound);
    const csMHD = mhdSoundSpeed(TeEv, TiEv, Z, mu, gammaElectron, gammaIon);
    const vTe   = electronThermalSpeed(TeEv);
    const vTi   = ionThermalSpeed(TiEv, mu);

    const betaI = speciesBeta(ni, TiEv, B);
    const betaE = speciesBeta(ne, TeEv, B);
    const betaTotal = betaI + betaE;

    const regime = kawRegime(betaE, mu);
    const kawRatio = regime.ratio;
    const kawState = regime.label;

    // rho_s = C_s^IA(gamma=1) / Omega_ci
    const rhos = ionSoundSpeed(TeEv, Z, mu, 1) / wci;

    // Taylor-shifted spectral break at the ion gyroscale (sub-ion break)
    const Vsw = VswKms * 1e3;
    const fbrk = Vsw > 0 ? Vsw / (2 * Math.PI * ri) : 0;

    return {
      fce, fci, fpe, fpi, fLH,
      lDe, ri, re, di, de,
      vA, csIA, csMHD, vTe, vTi,
      betaI, betaE, betaTotal,
      kawRatio, kawState,
      rhos, fbrk,
    };
  }

  global.PlasmaPhysics = {
    constants,
    conversions,
    electronGyroAngular,
    ionGyroAngular,
    electronPlasmaAngular,
    ionPlasmaAngular,
    lowerHybridAngular,
    electronDebyeLength,
    ionGyroradius,
    electronGyroradius,
    ionInertialLength,
    electronInertialLength,
    alfvenSpeed,
    ionSoundSpeed,
    mhdSoundSpeed,
    electronThermalSpeed,
    ionThermalSpeed,
    speciesBeta,
    electronIonMassRatio,
    kawRegime,
    alfvenEquivalentVelocity,
    plasmaState,
  };
})(window);
