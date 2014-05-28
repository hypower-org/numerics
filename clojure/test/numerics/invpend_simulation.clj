(ns numerics.invpend-simulation
  (:require [numerics.odesolve :as ode]
            [incanter.core :as inc]
            [incanter.charts :as charts]
            [incanter.pdf :as pdf]
            [clojure.core.matrix :as mat]
))

(mat/set-current-implementation :vectorz)

; This file demonstrates the simulation of a linearized inverted pendulum model using Euler.
; The model comes from Modern Control Systems, by Dorf and Bishop, 12th edition, pp. 186-187.
(def m-ball 0.2) ; mass of ball kg)
(def m-cart 1.0) ; mass of cart kg)
(def g 9.8) ; gravity
(def l 0.3) ; length of the rod (m)

(defn inv-pendulum-model
  [x t]
  (let [K (vector 1 1 1 1)
        x1 (first x)
        x2 (second x)
        x3 (nth x 2)
        x4 (nth x 3)
        u (mat/mul -1 (mat/dot K x))]
    (vector
      x2 ; first state equation
      (- (* (/ 1 m-cart) u) (* (/ (* m-ball g) m-cart) x3))  ; second state equation
      x4; third state equation
      (- (* (/ g l) x3) (* (/ 1 (* m-cart l)) u) ); fourth state equation
       )
    ))