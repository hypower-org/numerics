(ns numerics.odesolve)

(require '[clojure.core.matrix :as mat])
(mat/set-current-implementation :vectorz)

(defn my-round
  "A rounding function that rounds n to its closest value of dt."
  [n dt]
  (* (Math/round (/ n dt)) dt)
  )

(defn forward-euler
  "A function that numerically integrates a first order diffeq of the form xdot=f(x,t) forward in time over the 
interval t0 to tf."
  [diffeq x0 t0 tf delT]
  (loop [x x0 t t0 xs [] ts []]
    (if (> t tf)
      {:state xs :time ts}
      (recur (mat/add x (mat/mul (diffeq x t) delT)) (+ t delT) (conj xs x) (conj ts t)))))

(defn euler
  "A function that numerically integrates a first order diffeq of the form xdot=f(x,t) in time over the 
interval t0 to tf. If t0 > tf, then this solver operates in reverse."
  [diffeq x0 t0 tf delT]
  (let [op (if (> t0 tf)
             mat/sub
             mat/add)
        comp (if (> t0 tf)
             <=
             >=)]
    (loop [x x0 t t0 xs [] ts [] n 0]
      ;(println "num of xs: " (count xs))
      ;(println "time step: " (my-round t delT))
      (if (comp (my-round t delT) tf)
        {:state xs :time ts}
        ;(do 
          ;(println "recurring again..." n)
          (recur (op x (mat/mul (diffeq x t) delT)) (op t delT) (conj xs x) (conj ts t) (inc n) )
        ;  )
        )
      )
    )
  )

(defn forward-rk4
  "A function that applies 4th order runge kutta to numerically integrate a first order differential equation of 
the form xdot=f(x,t) over the interval t0 to tf."
  [diffeq x0 t0 tf delT]
  (loop [x x0 t t0 xs [] ts []]
    (if (> t tf)
      {:state xs :time ts}
      (let [h2 (* 0.5 delT)
            k1 (diffeq x t)
            k2 (diffeq (mat/add x (mat/mul k1 h2)) (+ t h2))
            k3 (diffeq (mat/add x (mat/mul k2 h2)) (+ t h2)) 
            k4 (diffeq (mat/add x (mat/mul k3 delT)) (+ t delT) ) ]
        (recur (mat/add x (mat/mul (mat/add k1 (mat/mul k2 2) (mat/mul k3 2) k4) (/ delT 6)) ) (+ t delT) (conj xs x) (conj ts t) ))
      )))
