(ns numerics.integration)

(defn forward-euler
  "A function that numerically integrates a diffeq (a scalar differential equation of the form xdot=f(x,t))
forward in time over the interval t0 to tf."
  [diffeq x0 t0 tf delT]
  (loop [x x0 t t0 xs [] ts []]
    (if (> t tf)
      {:state xs :times ts}
      (recur (+ x (* delT (diffeq x t))) (+ t delT) (conj xs x) (conj ts t)))))

(defn forward-rk4
  "A function that applies 4th order runge kutta to numerically integrate a scalar differential equation of 
the form xdot=f(x,t) over the interval t0 to tf."
  [diffeq x0 t0 tf delT]
  (loop [x x0 t t0 xs [] ts []]
    (if (> t tf)
      {:state xs :times ts}
      (let [h2 (* 0.5 delT)
            k1 (diffeq x t)
            k2 (diffeq (+ x (* h2 k1)) (+ t h2))
            k3 (diffeq (+ x (* h2 k2)) (+ t h2)) 
            k4 (diffeq (+ x (* delT k3)) (+ t delT) ) ]
        (recur (+ x (* (/ delT 6) (+ k1 (* 2 k2) (* 2 k3) k4))) (+ t delT) (conj xs x) (conj ts t) ))
      )))
