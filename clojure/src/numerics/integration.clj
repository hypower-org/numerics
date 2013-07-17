(ns numerics.integration)

; Note: this code came from the blog entry http://mattiasholmqvist.se/2009/10/integral-calculation-in-clojure/
; It's implementation is more "functional" than I was originally thinking.

(defn sub-series [f a h n index-filter]
  (map #(f (+ a (* h %))) (filter index-filter (range 1 n))))

(defn simpson [f a b n]
  (let [h (/ (- b a ) n)]
  (letfn [
      (sub-series [index-filter]
        (map #(f (+ a (* h %))) (filter index-filter (range 1 n))))]
    (* (/ h 3) (+
    (reduce + (map #(* % 2) (sub-series even?)))
    (reduce + (map #(* % 4) (sub-series odd?)))
    (f a)
    (f b))))))
  