use num_traits::{
    identities::{Zero, One},
    cast::FromPrimitive,
};

// use ordered_float::{Float, OrderedFloat};
/// Concrete computer representation of a real number 
type Real = f64; //OrderedFloat<f64>;

/// Size of the blind strategy
const M: usize = 30;

/// Detailed performance of the blind strategy. 
///
/// The guarantee is at least the minimum of the return vector.
fn detailed_performance(alpha: &[Real; M]) -> [Real; M + 1] {
    assert_eq!(alpha.len(), M);

    // Initialize
    let m = Real::from_usize(M).expect("M is real");
    let m_inv = m.recip();
    let one = Real::one();
    let zero = Real::zero();
    
    // Detailed factor 
    // g_{m,p}(k), where p = alpha_1
    let g = {
        let mut g = [zero; M];
        let p = alpha[0];
        for k in 1..=M {
            if k <= M/2 {
                g[k-1] = ( one - Real::from_usize(k).unwrap() * m_inv * ( one - p ) ).recip();
            }
            else {
                g[k-1] = Real::from_usize(2).unwrap() / ( one + p );
            }
        }
        g
    };

    // Computing garantee
    // f_j( \alpha_1, ..., \alpha_m )
    
    // Initialize
    let mut f = [zero; M + 1];
    // j = 1
    f[0] = (1..=M).map(|k| {
        alpha[0..(k-1)].iter().product::<Real>().powf(m_inv) 
            * (one - alpha[k-1].powf(m_inv))
            / -alpha[k-1].ln()
    }).sum();

    // j = m+1
    f[M] = one - alpha.iter().sum::<Real>() * m_inv;


    // j in {2, ..., m}
    // Initial building blocks
    let alpha_cumsum: [Real; M + 1] = {
        let mut acc = [zero; M + 1];
        for (i, a) in alpha.iter().enumerate() {
            acc[i + 1] = acc[i] + a;
        }
        acc
    };
    let alpha_cumprod = {
        let mut acc = [one; M + 1];
        for (i, a) in alpha.iter().enumerate() {
            acc[i + 1] = acc[i] * a;
        }
        acc
    };
    // Define for each j
    for j in 2..=M {
        let aux = (j..=M).map(|k| {
            let a_k = alpha[k-1];
            alpha_cumprod[k-1].powf(m_inv) * g[k-2] * (one - a_k.powf(m_inv)) / -a_k.ln()
        });
        f[j-1] = (Real::from_usize(j-1).expect("small usize is real") - alpha_cumsum[j-1]) * m_inv / (one - alpha[j-1])
            + aux.into_iter().sum::<Real>()
            ;
    }

    f
}

fn is_admissible(alpha: &[Real; M]) -> bool {
    // domain
    alpha.iter().all(|&a| Real::zero() <= a && a <= Real::one()) 
    &&
    // decreasing
    alpha.windows(2).all(|window: &[Real]| window[0] >= window[1])
} 

/// Problem cost function
fn cost(x: &[Real; M]) -> Real {
    if is_admissible(x) {
        - detailed_performance(x).into_iter().reduce(|a, b| a.min(b)).expect("There are elements")
    } else {
        Real::zero() // an undesirable value
    }
}

// Copyright 2018-2022 argmin developers
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

use argmin::core::observers::{ObserverMode, SlogLogger};
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::simulatedannealing::{Anneal, SATempFunc, SimulatedAnnealing};
use rand::distributions::Uniform;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::sync::{Arc, Mutex};

struct Problem {
    /// Random number generator. We use a `Arc<Mutex<_>>` here because `ArgminOperator` requires
    /// `self` to be passed as an immutable reference. This gives us thread safe interior
    /// mutability.
    rng: Arc<Mutex<Xoshiro256PlusPlus>>,
}

impl Problem {
    /// Constructor
    pub fn new() -> Self {
        Self {
            rng: Arc::new(Mutex::new(Xoshiro256PlusPlus::from_entropy())),
        }
    }
}

impl CostFunction for Problem {
    type Param = [Real; M];
    type Output = Real;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        Ok(cost(param))
    }
}

impl Anneal for Problem {
    type Param = [Real; M];
    type Output = [Real; M];
    type Float = Real;

    /// Anneal a parameter vector
    fn anneal(&self, param: &Self::Param, temp: Self::Float) -> Result<Self::Output, Error> {
        let mut param_n = param.clone();
        let mut rng = self.rng.lock().unwrap();
        let distr = Uniform::from(0..M);
        // Perform modifications to a degree proportional to the current temperature `temp`.
        for _ in 0..(temp.floor() as u64 + 1) {
            // Compute random index of the parameter vector using the supplied random number
            // generator.
            let idx = rng.sample(distr);

            // Compute random number in [-0.1, 0.1].
            let val = rng.sample(Uniform::new_inclusive(-0.1, 0.1));

            // modify previous parameter value at random position `idx` by `val`
            param_n[idx] += val;

            // check if bounds are violated. If yes, project onto bound.
            let lower_bound = {
                if idx == M - 1 {
                    (Real::zero() + param[idx]) * 0.5
                } else {
                    (param[idx + 1] + param[idx]) * 0.5
                }
            };
            let upper_bound = {
                if idx == 0 {
                    (Real::one() + param[idx]) * 0.5
                } else {
                    (param[idx - 1] + param[idx]) * 0.5
                }
            };
            param_n[idx] = param_n[idx].clamp(lower_bound, upper_bound);
        }
        Ok(param_n)
    }
}

fn run() -> Result<(), Error> {
    // Define cost function
    let operator = Problem::new();

    // Define initial parameter vector
    let init_param = [0.6270308410153931, 0.5684440914113474, 0.559466573817109, 0.5526891447996819, 0.5231248689913566, 0.518153567023222, 0.5141411139996803, 0.4768197283668264, 0.4435067420432881, 0.40337260112455486, 0.3993815651459096, 0.39626071800163687, 0.36552486861115796, 0.32662276834260257, 0.3080500632384776, 0.30520941038458327, 0.28791513316511735, 0.2876326328118053, 0.27511205040168396, 0.2689725713667366, 0.22908872176151002, 0.22829448594342994, 0.19833841963213145, 0.19006775086081584, 0.16161612625716001, 0.15301096038597622, 0.12937549169356408, 0.09210356754762833, 0.07471144203823056, 0.03364559149073742];
    // [0.5; M];

    // Define initial temperature
    let temp = 15.0;

    // Set up simulated annealing solver
    // An alternative random number generator (RNG) can be provided to `new_with_rng`:
    // SimulatedAnnealing::new_with_rng(temp, Xoshiro256PlusPlus::from_entropy())?
    let solver = SimulatedAnnealing::new(temp)?
        // Optional: Define temperature function (defaults to `SATempFunc::TemperatureFast`)
        .with_temp_func(SATempFunc::Boltzmann)
        /////////////////////////
        // Stopping criteria   //
        /////////////////////////
        // Optional: stop if there was no new best solution after 1000 iterations
        // .with_stall_best(1000)
        // Optional: stop if there was no accepted solution after 1000 iterations
        // .with_stall_accepted(1000)
        /////////////////////////
        // Reannealing         //
        /////////////////////////
        // Optional: Reanneal after 1000 iterations (resets temperature to initial temperature)
        // .with_reannealing_fixed(1000)
        // Optional: Reanneal after no accepted solution has been found for `iter` iterations
        // .with_reannealing_accepted(500)
        // Optional: Start reannealing after no new best solution has been found for 800 iterations
        // .with_reannealing_best(800)
        ;

    /////////////////////////
    // Run solver          //
    /////////////////////////
    let res = Executor::new(operator, solver)
        .configure(|state| {
            state
                .param(init_param)
                // Optional: Set maximum number of iterations (defaults to `std::u64::MAX`)
                .max_iters(1_000_000)
                // Optional: Set target cost function value (defaults to `std::f64::NEG_INFINITY`)
                .target_cost(-0.66975)
        })
        // Optional: Attach a observer
        .add_observer(SlogLogger::term(), ObserverMode::NewBest)
        .run()?;

    // Wait a second (lets the logger flush everything before printing again)
    std::thread::sleep(std::time::Duration::from_secs(1));

    // Print result
    println!("{res}");
    Ok(())
}

fn main() {
    if let Err(ref e) = run() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

// Previous results
//
// OptimizationResult:
//     Solver:        Simulated Annealing
//     param (best):  [0.6278188845507426, 0.618119857559708, 0.5775397759864143, 0.5760824650846582, 0.5659076677607481, 0.5391280119862243, 0.4686533337721845, 0.4380534045493786, 0.38270363958040343, 0.3780017030619437, 0.36348954705469966, 0.36299697081961746, 0.3581817109315925, 0.3523353070031661, 0.3218866288990676, 0.3132217673931533, 0.30877071449951743, 0.30289440938600287, 0.2910407407230008, 0.2642401436622602, 0.261491284085258, 0.24023833726087257, 0.20842835823235834, 0.19671917915790038, 0.16171196167609014, 0.15128957424863582, 0.1414531940538472, 0.09313092843162338, 0.05857758909259969, 0.04721290886241332]
//     cost (best):   -0.6674019705926073
//     iters (best):  5945
//     iters (total): 10000
//     termination:   Maximum number of iterations reached
//     time:          2.0372896s
//
// OptimizationResult:
//     Solver:        Simulated Annealing
//     param (best):  [0.649118189329632, 0.6258560004138201, 0.5813845104245882, 0.5587412098676289, 0.5090695172569351, 0.499566547804326, 0.47211024815010505, 0.4317808700503536, 0.422671360401653, 0.407369580377677, 0.3989162353145009, 0.37440338847703786, 0.37252157062582403, 0.3672163219302236, 0.3205056733594432, 0.3008392122568613, 0.2973644876549921, 0.28377987006804267, 0.27874632008740896, 0.2635507653768294, 0.2508656421322833, 0.24005671162121492, 0.21371263776087523, 0.1877267533713515, 0.16132365614563113, 0.1496388011609852, 0.13703729538482473, 0.09959643195006081, 0.0806388808317893, 0.031037362999340568]
//     cost (best):   -0.6677617982471253
//     iters (best):  7023
//     iters (total): 100000
//     termination:   Maximum number of iterations reached
//     time:          5.8481163s
//
// OptimizationResult:
//     Solver:        Simulated Annealing
//     param (best):  [0.6275888862613779, 0.6192199516101676, 0.6139885808216932, 0.5847786128356477, 0.5255391965662364, 0.5032224561732854, 0.4701083149708849, 0.4603667239656413, 0.4357550526242747, 0.4227220092955919, 0.4142186466090482, 0.3835792173897641, 0.372990274065627, 0.3194576855809654, 0.3101466609453093, 0.29965249828483975, 0.2890032071765063, 0.26760609683217934, 0.26518443371263356, 0.25275541884564073, 0.23453841033617273, 0.22111862025569579, 0.2127536660565621, 0.19238304587796332, 0.16733908700390498, 0.16190980556978107, 0.1241869192160632, 0.09050797035246794, 0.06025525974859788, 0.031084577445222337]
//     cost (best):   -0.6682451292282547
//     iters (best):  40546
//     iters (total): 1000000
//     termination:   Maximum number of iterations reached
//     time:          80.4118751s
//
// OptimizationResult:
//     Solver:        Simulated Annealing
//     param (best):  [0.6270308410153931, 0.6184440914113474, 0.559466573817109, 0.5526891447996819, 0.5231248689913566, 0.518153567023222, 0.5141411139996803, 0.4768197283668264, 0.4435067420432881, 0.40337260112455486, 0.3993815651459096, 0.39626071800163687, 0.36552486861115796, 0.32662276834260257, 0.3080500632384776, 0.30520941038458327, 0.28791513316511735, 0.2876326328118053, 0.27511205040168396, 0.2689725713667366, 0.22908872176151002, 0.22829448594342994, 0.19833841963213145, 0.19006775086081584, 0.16161612625716001, 0.15301096038597622, 0.12937549169356408, 0.09210356754762833, 0.07471144203823056, 0.03364559149073742]
//     cost (best):   -0.6683340031470274
//     iters (best):  932905
//     iters (total): 1000000
//     termination:   Maximum number of iterations reached
//     time:          74.7414323s
// 
// OptimizationResult:
//     Solver:        Simulated Annealing
//     param (best):  [0.6488377139147816, 0.6286748255066426, 0.5959259780730781, 0.581146287337283, 0.5127584720607117, 0.49709499001209534, 0.4800210595809111, 0.46880034959807604, 0.4383713763928919, 0.4073292082536808, 0.36363890895975204, 0.35956533568360116, 0.34557391150818206, 0.34189193603737095, 0.33733738779663586, 0.3292091977620634, 0.3264539322982454, 0.2821942725093516, 0.263229188864773, 0.24683708037169247, 0.24287176137061234, 0.21895174321906052, 0.21534660761148774, 0.1877662629954155, 0.1580734012651116, 0.14218799779672164, 0.11287582062790322, 0.10297556903027627, 0.07025581930535822, 0.04094300669695767]
//     cost (best):   -0.6684144764570266
//     iters (best):  268808
//     iters (total): 1000000
//     termination:   Maximum number of iterations reached
//     time:          76.6206991s