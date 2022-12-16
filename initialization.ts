export function initialization() {
    const inletPressure = 25000000; // Па
    const inletTemperature = 300; // K
    const atmospherePressure = 100000; // Па

    return {
        boundaryConditions: {
            P1: inletPressure,
            Pa: atmospherePressure,
            T1: inletTemperature,
            Ta: inletTemperature,
            P2: inletPressure,
            Xp: 0,
            _Xp: 0,
            T2: inletTemperature,
            Vp: 0,
            _Vp: 0,
            P3: inletPressure,
            T3: inletTemperature,
            P4: atmospherePressure,
            T4: inletTemperature,
            P5: atmospherePressure,
            T5: inletTemperature,
        },

        T_init: inletTemperature,
        pvh: inletPressure,
        pa: atmospherePressure,
        ph: atmospherePressure,
        Ta: inletTemperature,
        t: 0,
        w: 0,
        topen: 0,
        tclose: 0,
        trab: 0,
        m2: 0.7,
        m4: 0,
        u: 0,
    } as const;
}
