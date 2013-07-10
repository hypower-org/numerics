(ns numerics.odesolve)

(require '[mikera.vectorz.core :as vz])

(defn forward-euler
  "A function that numerically integrates a first order diffeq of the form xdot=f(x,t) forward in time over the 
interval t0 to tf."
  [diffeq x0 t0 tf delT]
  (let [add (if (vz/vec? x0) ; check if the initial condition is a vectorz vector, or just a scalar!
              vz/+
              +)
        mult (if (vz/vec? x0) ; check if the initial condition is a vectorz vector, or just a scalar!
              vz/scale
              *)]
  (loop [x x0 t t0 xs [] ts []]
    (if (> t tf)
      {:state xs :time ts}
      (recur (add x (mult (diffeq x t) delT)) (+ t delT) (conj xs x) (conj ts t))))))

(defn forward-rk4
  "A function that applies 4th order runge kutta to numerically integrate a first order differential equation of 
the form xdot=f(x,t) over the interval t0 to tf."
  [diffeq x0 t0 tf delT]
  (let [add (if (vz/vec? x0) ; check if the initial condition is a vectorz vector, or just a scalar!
              vz/+
              +)
        mult (if (vz/vec? x0) ; check if the initial condition is a vectorz vector, or just a scalar!
              vz/scale
              *)]
  (loop [x x0 t t0 xs [] ts []]
    (if (> t tf)
      {:state xs :time ts}
      (let [h2 (* 0.5 delT)
            k1 (diffeq x t)
            k2 (diffeq (add x (mult k1 h2)) (+ t h2))
            k3 (diffeq (add x (mult k2 h2)) (+ t h2)) 
            k4 (diffeq (add x (mult k3 delT)) (+ t delT) ) ]
        (recur (add x (mult (add k1 (mult k2 2) (mult k3 2) k4) (/ delT 6)) ) (+ t delT) (conj xs x) (conj ts t) ))
      ))))
