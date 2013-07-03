(ns numerics.integration-test
  (:use clojure.test
        numerics.integration))

(require '[mikera.vectorz.core :as vz])
(use '(incanter core charts))

; Here we simulate the vector differential equation dot{x}_1 = -x_1 + x_2; dot{x}_2 = -x_2:
(defn vector-ode [x t]
  (vz/vec [(+ (* -1 (first x)) (second x)) (* -1 (second x)) ]))

(let [x0 (vz/vec [2 -1]) ; set initial condition to some vector
      results (forward-euler vector-ode x0 0 1 0.01) ; solve numerically using Euler's method
      x1 (map first (:state results)) ; extract the states
      x2 (map second (:state results))
      result-plot (xy-plot (:time results) x1)] ; now use Incanter to display our wonderful results...
  (do
    (add-lines result-plot (:time results) x2)
    (view result-plot)))
  