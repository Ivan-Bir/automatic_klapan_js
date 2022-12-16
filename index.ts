const fs = require('fs'); // модуль для работы с файловой системой
import { initialization } from "./initialization";
import { getConstants } from "./getConstants";
import {MathModelStateCoeffs, MathModelCharacteristics, MathModelObject, RK_factors, Solver, State} from "./types";

const outTableFileName = 'outTable.txt';
fs.writeFileSync(outTableFileName, '');

const props = getConstants();
const inits = initialization();

let t = inits.t;

let currParameters: MathModelObject = {
    functionalParams: {...inits.boundaryConditions}, // копирование объекта инициализации для 0-й итерации
    options: { // Параметры, необходимые уже на первой интерации расчёта
        m2: inits.m2,
        m4: inits.m4,
        f1: props.m1*areaCirc(props.d1),
        f3: props.m3*areaCirc(props.d3),
        f6: props.m6*areaCirc(props.d6),
        f8: props.m8*areaCirc(props.d8),
        fkp: areaCirc(0.0175),
    }
};

let state: State = { // Переменная состояния рабочего процесса
    isClosed: 1,
    isOpened: 0,
    isRabProcess: 0,
    fromOpenToClose: 0,
    fromCloseToOpen: 0,
    isDone: 0,
};

let j = 0;
while(!state.isOpened) {
    const RK_0 = calcLeftBoundaryRungeKutta(currParameters, MathModelOpenIteration);
    const RK_1 = calcFirstCenterRungeKutta(currParameters, props.dt, RK_0, MathModelOpenIteration);
    const RK_2 = calcSecondaryCenterRungeKutta(currParameters, props.dt, RK_1, MathModelOpenIteration);
    const RK_3 = calcRightBoundaryRungeKutta(currParameters, props.dt, RK_2, MathModelOpenIteration);

    currParameters.functionalParams = updateParameters(currParameters.functionalParams, RK_0, RK_1, RK_2, RK_3, props.dt);
    if (state.isClosed) {
        currParameters.options.m2 = 0;
        currParameters.options.m4 = 0.7;
    }

    state = updateState(state, currParameters.functionalParams);

    if (state.isOpened) {
        currParameters.options.m2 = 0.7;
        currParameters.options.m4 = 0;
    }

    if (currParameters.functionalParams.Xp <= 0 && currParameters.functionalParams.Vp <= 0) {
        currParameters.functionalParams.Vp = 0;
        currParameters.functionalParams._Xp = props.h;
        currParameters.functionalParams._Vp = 0;
    }

    if (currParameters.functionalParams.Xp > props.h ) {
        currParameters.functionalParams.Vp = 0;
        currParameters.functionalParams._Xp = 0;
        currParameters.functionalParams._Vp = 0;
    }

    currParameters.functionalParams.Xp = Math.max(0,Math.min(currParameters.functionalParams.Xp, props.h));
    currParameters.functionalParams._Xp = Math.max(0,Math.min(currParameters.functionalParams._Xp, props.h));

    const {P1,T1,P2,T2,P3,T3,P4,T4,P5,T5,Pa,Ta,Xp,Vp,_Xp,_Vp} = currParameters.functionalParams;

    if (j % 10 === 0) { // Каждые 10 итераций вывод в файл
        fs.appendFileSync(outTableFileName, `${t.toFixed(5)}  ${P1.toFixed(0)}  ${P2.toFixed(0)}  ${P3.toFixed(0)}  ${P4.toFixed(0)}  ${P5.toFixed(0)}  ${T1.toFixed(3)}  ${T2.toFixed(3)}  ${T3.toFixed(3)}  ${T4.toFixed(3)}  ${T5.toFixed(3)}  ${Pa}  ${Ta}  ${Xp.toFixed(6)}  ${Vp.toFixed(6)}  \n`);

        if (j % 100 === 0 ) { //Каждые 100 итераций печать на экран
            console.clear();
            console.log(`t=${t.toFixed(4)}\nP1=${P1.toFixed(0)}\tP2=${P2.toFixed(0)}\tP3=${P3.toFixed(0)}\tP4=${P4.toFixed(0)}\tP5=${P5.toFixed(0)}
T1=${T1.toFixed(3)}\tT2=${T2.toFixed(3)}\tT3=${T3.toFixed(3)}\tT4=${T4.toFixed(3)}\tT5=${T5.toFixed(3)}
Xp=${Xp.toFixed(5)}\tVp=${Vp.toFixed(5)}\t_Xp=${_Xp.toFixed(5)}\t_Vp=${_Vp.toFixed(5)}`);
            j = 0;
        }
    }
    t += props.dt;
    ++j;

    if (state.isDone) {
        break;
    }
};

console.log('SUCCESS');

// Вычисляет массовый расход
// Проходная площадь f [м^2]
// Давление на входе p1 [Па]
// Давление на выходе p2 [Па]
// Температура в области T [К]
function MassFlow(f:number, p1:number, p2:number, T:number) {
    const { k, R, bk } = props;
    if (p2/p1 <= bk) {
        return f*p1*Math.sqrt(2*k/((k-1)*R*T)*Math.abs(Math.pow(bk,(2/k))-Math.pow(bk,((k+1)/k))));
    }

    return f*p1*Math.sqrt(2*k/((k-1)*R*T)*(Math.pow(bk,(2/k))-Math.pow(bk,((k+1)/k))));
}

// Вычисляет подогрев
// Tемпература T
// s ??
function dQ(T: number, s:number) {
    const { alpha } = props;
    const { T_init } = inits;
    return alpha * s * (T - T_init);
}

// Вычисляет площадь круга
function areaCirc(diameter: number) {
    return diameter*diameter*Math.PI/4;
}

// Определяет знак;
function sig(a:number, b:number):number {
    if (a >= b) {
        return 1;
    }

    return -1;
}

// Функция Математической модели.
// Принимает динамические характеристики предыдущего состояния модели и управляющие параметры модели
// Возвращает коэффициенты для формулы Рунге-Кутта
function MathModelOpenIteration(params: MathModelCharacteristics, mutableOptions: MathModelStateCoeffs): RK_factors {
    const {d2, d4} = props;
    mutableOptions.f2 = mutableOptions.m2*areaCirc(d2);
    mutableOptions.f4 = mutableOptions.m4*(areaCirc(d4)-areaCirc(0.00075)); //?

    const {pvh, pa, Ta} = inits;

    mutableOptions.pa_vh1 = Math.max(pvh, params.P1);
    mutableOptions.pb_vh1 = Math.min(pvh,params.P1);
    mutableOptions.pa_12 = Math.max(params.P1,params.P2);
    mutableOptions.pb_12 = Math.min(params.P1,params.P2);
    mutableOptions.pa_23 = Math.max(params.P2,params.P3);
    mutableOptions.pb_23 = Math.min(params.P2,params.P3);
    mutableOptions.pa_2a = Math.max(params.P2,pa);
    mutableOptions.pb_2a = Math.min(params.P2,pa);
    mutableOptions.pa_14 = Math.max(params.P1,params.P4);
    mutableOptions.pb_14 = Math.min(params.P1,params.P4);
    mutableOptions.pa_4v = Math.max(params.P4,pa);
    mutableOptions.pb_4v = Math.min(params.P4,pa);
    mutableOptions.pa_45 = Math.max(params.P4,params.P5);
    mutableOptions.pb_45 = Math.min(params.P4,params.P5);
    mutableOptions.pa_5a = Math.max(params.P5,pa);
    mutableOptions.pb_5a = Math.min(params.P5,pa);

    mutableOptions.Tvh1 = (pvh >= params.P1)? Ta : params.T1; // (условие)? <выполняется если True>: <выполняется в противном случае>
    mutableOptions.T12 = (params.P1 >= params.P2)? params.T1 : params.T2;
    mutableOptions.T23 = (params.P2 >= params.P3)? params.T2 : params.T3;
    mutableOptions.Ta2 = (params.P2 >= pa)? params.T2 : params.Ta;
    mutableOptions.T14 = (params.P1 >= params.P4)? params.T1 : params.T4;
    // mutableOptions.Tv4 = (pa >= params.P4)? params.Ta : params.T4;
    mutableOptions.Tv4 = (params.P4 >= pa)? params.T4 : params.Ta;
    mutableOptions.T45 = (params.P4 >= params.P5)? params.T4 : params.T5;
    mutableOptions.Ta5 = (pa >= params.P5)? Ta : params.T5;

    if (params.Xp < 0 && params.Vp < 0) {
        params.Vp = 0;
        params._Xp = props.h;
        params._Vp = 0;
    }

    if (params.Xp > props.h ) {
        params.Vp = 0;
        params._Xp = 0;
        params._Vp = 0;
    }

    params.Xp = Math.max(0, Math.min(params.Xp, props.h));
    params._Xp = Math.max(0, Math.min(params._Xp, props.h));

    mutableOptions.fkl = (params.Xp <= 0)? mutableOptions.fkp - areaCirc(0.0095): mutableOptions.fkp - areaCirc(0.004);
    mutableOptions.gvh1 = sig(pvh, params.P1)*MassFlow(mutableOptions.f1,
                                                        mutableOptions.pa_vh1,
                                                        mutableOptions.pb_vh1,
                                                        mutableOptions.Tvh1);

    mutableOptions.g1vh = -mutableOptions.gvh1;

    mutableOptions.g12 = sig(params.P1, params.P2)*MassFlow(mutableOptions.f2,
                                                            mutableOptions.pa_12,
                                                            mutableOptions.pb_12,
                                                            mutableOptions.T12);

    mutableOptions.g21 = mutableOptions.g12;

    mutableOptions.g23 = sig(params.P2, params.P3)*MassFlow(mutableOptions.f3,
                                                            mutableOptions.pa_23,
                                                            mutableOptions.pb_23,
                                                            mutableOptions.T23);

    mutableOptions.g32 = -mutableOptions.g23;

    const {m7, m5, d5, d7, h } = props;
    mutableOptions.f5 = m5*Math.PI*d5*params.Xp;
    mutableOptions.f7 = m7*Math.PI*d7*(h-params.Xp);

    mutableOptions.g2a = sig(params.P2,pa)*MassFlow(mutableOptions.f4,
                                                    mutableOptions.pa_2a,
                                                    mutableOptions.pb_2a,
                                                    mutableOptions.Ta2);

    mutableOptions.ga2 = -mutableOptions.g2a;

    mutableOptions.g14 = sig(params.P1, params.P4)*MassFlow(mutableOptions.f5,
                                                            mutableOptions.pa_14,
                                                            mutableOptions.pb_14,
                                                            mutableOptions.T14);

    mutableOptions.g41 = -mutableOptions.g14;

    mutableOptions.g4v = sig(params.P4, pa)*MassFlow(mutableOptions.f6,
                                                    mutableOptions.pa_4v,
                                                    mutableOptions.pb_4v,
                                                    mutableOptions.Tv4);

    mutableOptions.gv4 = -mutableOptions.g4v;

    mutableOptions.g45 = sig(params.P4, params.P5)*MassFlow(mutableOptions.f7,
                                                            mutableOptions.pa_45,
                                                            mutableOptions.pb_45,
                                                            mutableOptions.T45);

    mutableOptions.g54 = -mutableOptions.g45;

    mutableOptions.g5a = sig(params.P5, pa)*MassFlow(mutableOptions.f8,
                                                    mutableOptions.pa_5a,
                                                    mutableOptions.pb_5a,
                                                    mutableOptions.Ta5);

    mutableOptions.ga5 = -mutableOptions.g5a;

    const {v30} = props;
    mutableOptions.v3 = v30-params.Xp*mutableOptions.fkp;
    mutableOptions.v3p = -params.Vp*mutableOptions.fkp;

    const {v10, s10, s30} = props;
    mutableOptions.v1 = v10+params.Xp*mutableOptions.fkl;
    mutableOptions.v1p = params.Vp*mutableOptions.fkl;

    mutableOptions.s1 = s10+params.Xp*Math.PI*0.0175; // ?
    mutableOptions.s3 = s30-params.Xp*Math.PI*0.0175; // ?

    const {k,R, s2, s4, s5, v2, v4, v5} = props;

    const {gvh1, g21, g41, v1p, s1, v1, Tvh1, T12, T14} = mutableOptions;
    const CoefP1 = k*(R*Tvh1*gvh1+R*T12*g21+R*T14*g41-params.P1*v1p-(k-1)*dQ(params.T1,s1))/v1;
    const CoefT1 = params.T1/(params.P1*v1)*(v1*CoefP1+params.P1*v1p-R*(Tvh1*gvh1+T12*g21+T14*g41)-(k-1)*dQ(params.T1,s1));

    const {g12, g32, ga2, T23, Ta2} = mutableOptions;
    const CoefP2 = k*(R*T12*g12+R*T23*g32+R*Ta2*ga2-(k-1)*dQ(params.T2, s2))/v2;
    const CoefT2 = params.T2/(params.P2*v2)*(v2*CoefP2-R*(T12*g12+T23*g32+Tvh1*ga2)-(k-1)*dQ(params.T2,s2));

    const {g23, v3p, s3, v3} = mutableOptions;
    const CoefP3 = k*(R*T23*g23-params.P3*v3p-(k-1)*dQ(params.T3,s3))/v3;
    const CoefT3 = params.T3/(params.P3*v3)*(v3*CoefP3+params.P3*v3p-R*T23*g23-(k-1)*dQ(params.T3,s3));

    const {g14, gv4, g54,Tv4, T45} = mutableOptions;
    const CoefP4 = k*(R*T14*g14+R*Tv4*gv4+R*T45*g54-(k-1)*dQ(params.T4,s4))/v4;
    const CoefT4 = params.T4/(params.P4*v4)*(v4*CoefP4-R*(T14*g14+Tv4*gv4+T45*g54)-(k-1)*dQ(params.T4,s4));

    const {g45, ga5,Ta5} = mutableOptions;
    const CoefP5 = k*(R*T45*g45+R*Ta5*ga5-(k-1)*dQ(params.T5,s5))/v5;
    const CoefT5 = params.T5/(params.P5*v5)*(v5*CoefP5-R*(T45*g45+Ta5*ga5)-(k-1)*dQ(params.T5,s5));

    const CoefPa = 0;
    const CoefTa = 0;

    const {Cp, ktr, x0, M} = props;
    const CoefXp = params.Vp;
    const CoefVp = (mutableOptions.fkl*params.P1-mutableOptions.fkp*params.P3-ktr*params.Vp-Cp*(x0+params.Xp)
         + pa*areaCirc(0.016)-params.P4*(areaCirc(0.016)-areaCirc(0.004)))/M;
    const Coef_Xp = params._Vp;
    const Coef_Vp = -CoefVp;

    const updatedParameters: RK_factors = {
        P1: CoefP1,
        T1: CoefT1,
        P2: CoefP2,
        T2: CoefT2,
        P3: CoefP3,
        T3: CoefT3,
        P4: CoefP4,
        T4: CoefT4,
        Pa: CoefPa,
        Ta: CoefTa,
        P5: CoefP5,
        T5: CoefT5,
        Xp: CoefXp,
        Vp: CoefVp,
        _Xp: Coef_Xp,
        _Vp: Coef_Vp,
    };

    return updatedParameters;
}

function calcLeftBoundaryRungeKutta(currParameters: MathModelObject, solveIteration: Solver) : RK_factors {
    return solveIteration(currParameters.functionalParams, currParameters.options);
};

function calcFirstCenterRungeKutta(currParameters: MathModelObject, dt:number, C0: RK_factors, solveIteration: Solver) : RK_factors {
    const modifParams = {};
    for (const key in currParameters.functionalParams) {
        // увеличиваем значения каждой характеристики на C0*h/2
        (modifParams as any)[key] = (currParameters.functionalParams as any)[key] + (C0 as any)[key]*dt/2;
    }

    return solveIteration(modifParams as MathModelCharacteristics, currParameters.options);
};

function calcSecondaryCenterRungeKutta(currParameters: MathModelObject, dt:number, C1: RK_factors, solveIteration: Solver) : RK_factors {
    const modifParams = {};
    for (const key in currParameters.functionalParams) {
        // увеличиваем значения каждой характеристики на C1*h/2
        (modifParams as any)[key] = (currParameters.functionalParams as any)[key] + (C1 as any)[key]*dt/2;
    }

    return solveIteration(modifParams as MathModelCharacteristics, currParameters.options);
};

function calcRightBoundaryRungeKutta(currParameters: MathModelObject, dt:number, C2: RK_factors, solveIteration: Solver) : RK_factors {
    const modifParams = {};
    for (const key in currParameters.functionalParams) {
        // увеличиваем значения каждой характеристики на C2*h
        (modifParams as any)[key] = (currParameters.functionalParams as any)[key] + (C2 as any)[key]*dt;
    }

    return solveIteration(modifParams as MathModelCharacteristics, currParameters.options);
};

function updateParameters(params: MathModelCharacteristics, RK_0: RK_factors, RK_1: RK_factors, RK_2: RK_factors, RK_3: RK_factors, dt:number): MathModelCharacteristics {
    const updatedParams = {};
    for (const key in params) {
        // увеличиваем значения каждой характеристики на величину приращения по формуле Рунге-Кутта
        (updatedParams as any)[key] = (params as any)[key] + dt*(
            (RK_0 as any)[key] + 2*(RK_1 as any)[key] + 2*(RK_2 as any)[key] + (RK_3 as any)[key]
        )/6;
    }

    return updatedParams as MathModelCharacteristics;
}

function updateState(st:State, params:MathModelCharacteristics):State {
    if (st.isClosed) {
        return {
            isClosed: 0,
            isOpened: 0,
            isRabProcess: 0,
            fromOpenToClose: 0,
            fromCloseToOpen: 1,
            isDone: 0,
        };
    }

    if (0.9*params.P3 <= inits.pa && st.fromCloseToOpen) {
        return {
            isClosed: 0,
            isOpened: 0,
            isRabProcess: 1,
            fromOpenToClose: 0,
            fromCloseToOpen: 0,
            isDone: 0,
        };
    }

    if (0.95*params.P3 <= inits.pa && st.isRabProcess) {
        return {
            isClosed: 0,
            isOpened: 0,
            isRabProcess: 0,
            fromOpenToClose: 0,
            fromCloseToOpen: 0,
            isDone: 1,
        };
    }

    if (params.P3 > 0.9*inits.pvh && st.isOpened) {
        return {
            isClosed: 0,
            isOpened: 0,
            isRabProcess: 0,
            fromOpenToClose: 0,
            fromCloseToOpen: 0,
            isDone: 1,
        };
    }

    return st;
}
