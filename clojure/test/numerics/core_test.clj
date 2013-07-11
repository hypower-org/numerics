(ns numerics.integration-test
  (:use clojure.test
        numerics.odesolve))

(use '(incanter core charts))

; Here we simulate the vector differential equation dot{x}_1 = -x_1 + x_2; dot{x}_2 = -x_2:
(defn vector-ode [x t]
  (vec [(+ (* -1 (first x)) (second x)) (* -1 (second x)) ]))

(defn scalar-ode [x t]
  (* -0.8 x))

; vector ode integration ( 2 variables)
(time (let [x0 (vec [2 -1]) ; set initial condition to some vector
      results (forward-euler vector-ode x0 0 10 0.01) ; solve numerically using Euler's method
      x1 (map first (:state results)) ; extract the states
      x2 (map second (:state results))
      result-plot (xy-plot (:time results) x1)] ; now use Incanter to display our wonderful results...
  (do
    (add-lines result-plot (:time results) x2)
    (view result-plot))))

(time (let [x0 (vec [3 -1]) ; set initial condition to some vector
      results (forward-rk4 vector-ode x0 0 10 0.01) ; solve numerically using Euler's method
      x1 (map first (:state results)) ; extract the states
      x2 (map second (:state results))
      result-plot (xy-plot (:time results) x1)] ; now use Incanter to display our wonderful results...
  (do
    (add-lines result-plot (:time results) x2)
    (view result-plot))))

; scalar ode integration
(defn scalar-ode [x t]
  (* -0.8 x))

(time (let [x0 2 ; set initial condition to some vector
            results (forward-euler scalar-ode x0 0 4 0.01) ; solve numerically using Euler's method
            x1 (:state results) ; extract the states
            result-plot (xy-plot (:time results) x1)] ; now use Incanter to display our wonderful results...
        (view result-plot)))

(time (let [x0 2 ; set initial condition to some vector
            results (forward-rk4 scalar-ode x0 0 4 0.01) ; solve numerically using Euler's method
            x1 (:state results) ; extract the states
            result-plot (xy-plot (:time results) x1)] ; now use Incanter to display our wonderful results...
        (view result-plot)))

; plot the graphs to compare accuracy 
(let [x0 2
      euler-results (forward-euler scalar-ode x0 0 4 0.01)
      rk4-results (forward-rk4 scalar-ode x0 0 4 0.01)
      result-plot (xy-plot (:time euler-results) (:state euler-results))]
  (do 
    (add-lines result-plot (:time rk4-results) (:state rk4-results))
    (view result-plot)))
      
