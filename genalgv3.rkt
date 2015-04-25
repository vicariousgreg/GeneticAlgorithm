#lang plai-typed
(print-only-errors true)
(require (typed-in racket
                   [abs : (number -> number)]
                   [build-list : (number (number -> number) -> (listof number))]
                   [displayln : ('a -> void)]
                   [exact->inexact : (number -> number)]
                   [gensym : (-> symbol)]
                   [shuffle : ((listof 'a) -> (listof 'a))]))

(require (rename-in (typed-in racket [random : (-> number)])
                    [random random])
         (rename-in (typed-in racket [random : (number -> number)])
                    [random random-range]))

;;;;;;;;;;;;;;;;;;;;;;;; OVERVIEW ;;;;;;;;;;;;;;;;;;;;;;;;;
; The following is an implementation of a simple genetic
; algorithm including the following features:
;  Tournament selection
;  Crossover & mutation
;  Elitism
;  Generation capping

; Primary algorithm parameters include:
;  Population size
;  Tournament size
;  Mutation rate
;  Elitism flag
;  Generation cap
;  Print generations
; Secondary algorithm parameters include:
;  Individual size
;  Algorithm solution
;  Solution score
;  Max difference
;  Gene size

; Secondary algorithm parameters are implementation
; dependent, as they are particular to the type of
; individuals in the population.

; To customize the algorithm, the secondary parameters
; must be configured, and the following functions must
; be re-implemented:
;  gen-tour
;  calc-tour-length
;  crossover
;  mutate-tour
;;;;;;;;;;;;;;;;;;;;;; V2 NOTES ;;;;;;;;;;;;;;;;;;;;;;;;;;
; Version 2 features integers for genes instead of bits.
; Fitness score is replaced with a difference score which
; is to be minimized by the algorithm.

; Additionally, a flag has been added for printing the
; generations and their best scores, and the main
; algorithm function has been improved to print the found
; solution and the number of generations it took to
; discover.
;;;;;;;;;;;;;;;;;;;;;; V3 NOTES ;;;;;;;;;;;;;;;;;;;;;;;;;;
; Version 3 features improvements for finding solutions
; to the travelling salesperson problem.  Graphs were
; added to represent TSP graphs of city connections and
; their weights.  Crossing over and mutation have been
; updated to operate on TSP tours.  Since an optimal
; solution is not provided, the generation cap is required
; for algorithm termination.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;; HELPERS ;;;;;;;;;;
; Retrieves the nth elementh of a list.
(define (nth [n : number] [l : (listof 'a)]) : 'a
  (cond
    [(= 0 n) (error 'nth "Cannot get 0th element of list")]
    [(empty? l) (error 'nth "List not long enough")]
    [else (if (= n 1) (first l) (nth (sub1 n) (rest l)))]))

(test (nth 5 (list 1 2 3 4 5)) 5)
(test (nth 1 (list 1 2 3)) 1)
(test (nth 3 (list 1 2 3)) 3)
(test/exn (nth 0 (list 1 2 3 4)) "nth")
(test/exn (nth 2 (list 1)) "nth")

; Consumes two lists of the same length and returns
;  a new list of tuples where each element of the new
;  tuple is a list containing both of the corresponding
;  elements from the original lists.
(define (zip [l : (listof 'a)]
             [r : (listof 'b)]) : (listof ('a * 'b))
  (cond
    [(and
      (or (empty? l) (empty? r))
      (not (and (empty? l) (empty? r))))
     (error 'zip "Lists must be the same length.")]
    [(empty? l) empty]
    [else (cons (values (first l) (first r))
      (zip (rest l) (rest r)))]))

(test (zip
       (list 1 2 3)
       (list 4 5 6))
      (list
       (values 1 4)
       (values 2 5)
       (values 3 6)))
(test (zip
       (list 1)
       (list 4))
      (list
       (values 1 4)))
(test (zip empty empty)
      empty)
(test/exn (zip (list 1) (list 1 2))
          "Lists must be the same length.")
(test/exn (zip (list 1) empty)
          "Lists must be the same length.")


;;;;;;;;;; PARAMETERS ;;;;;;;;;;

;;;;; Types ;;;;;
(define-type-alias City symbol)
(define-type-alias Tour (listof City))
(define-type-alias Population (listof Tour))

;;;;; Algorithm Parameters ;;;;;
; The size of the population
(define kPopulationSize 50)

; The number of cities in the graph
(define kGraphSize 100)

; The number of tours selected for a tournament
(define kTournamentSize 5)

; The rate of random mutation
(define kMutationRate 0.015)

; Generation cap
;  Caps the number of generations so that the best solution is
;  returned after a certain number of iterations.
; A generation cap of 0 means the algorithm is allowed to run
; until an acceptable solution is found.
(define kGenerationCap 1000)

; Print generations flag
;  Enables/disables printing of each generation's current best
(define kPrintGenerations false)


;;;;;;;;;; GRAPHS ;;;;;;;;;;
; Graph datatype
(define-type-alias Graph (hashof symbol (hashof symbol number)))

; City list generator.
; Generates a list of unique city names as symbols.
(define (generate-cities [size : number]) : (listof City)
  (cond
    [(= size 0) empty]
    [else (cons (gensym) (generate-cities (sub1 size)))]))

; Generates a list of cities and weights.
; Represents the weights of connections from the given city to each
; other city in the list.
(define (generate-connections [city : symbol] [cities : (listof City)])
  : (listof (symbol * number))
  (cond
    [(empty? cities) empty]
    [(symbol=? city (first cities))
     (generate-connections city (rest cities))]
    [else (cons (values (first cities) (random-range 20))
                (generate-connections city (rest cities)))]))

; Graph generator.
; Generates a graph of the given size and returns the list of cities
; and the graph of connections.
(define (generate-graph [size : number]) : ((listof City) * Graph)
  (let ([cities (generate-cities size)])
    (values cities 
            (hash (zip cities
                       (map (lambda (city)
                              (hash (generate-connections city cities)))
                            cities))))))

;;;;;;;;;; ALGORITHM MACHINERY ;;;;;;;;;;

;;;; 1. INITIALIZATION ;;;;;
; Generates a population randomly.
;   Population size : kPopulationSize
(define (gen-pop [cities : Tour]) : Population
  (gen-pop-helper kPopulationSize cities))

; Recursive helper for population generation.
;  Generates a population of 
(define (gen-pop-helper [size : number] [cities : Tour]) : Population
  (cond
    [(= size 0) empty]
    [else (cons (gen-tour cities)
                (gen-pop-helper (sub1 size) cities))]))

; Generates a tour.
(define (gen-tour [cities : Tour]) : Tour
  (shuffle cities))

;;;; 2. EVALUATION
; Evaluates the tour length of an tour from the graph.
(define (calc-tour-length [indiv : Tour] [graph : Graph]) : number
  (cond
    [(= 1 (length indiv)) 0]
    [else (let ([recur (calc-tour-length (rest indiv) graph)])
            (+ recur (some-v (hash-ref (some-v (hash-ref graph (first indiv)))
                                       (second indiv)))))]))

;;;; 3. SELECTION
; Gets the shortest tour of a population.
(define (get-best [pop : Population] [graph : Graph]) : Tour
  (best-helper (rest pop)
               (first pop)
               (calc-tour-length (first pop) graph) graph))

; Recursive helper for best tour.
(define (best-helper [pop : Population]
                     [best : Tour]
                     [best-difference : number]
                     [graph : Graph]) : Tour
  (cond
    [(empty? pop) best]
    [else (local ([define curr (first pop)]
                  [define curr-length (calc-tour-length curr graph)])
            (cond
              [(< curr-length best-difference)
               (best-helper (rest pop) curr curr-length graph)]
              [else (best-helper (rest pop) best best-difference graph)]))]))

; Performs a tournament and returns the winning tour.
;  Tournament size : kTournamentSize
(define (tournament [pop : Population] [graph : Graph]) : Tour
  (get-best (gen-tournament-pop pop kTournamentSize) graph))

; Generates a tournament population of the given size.
(define (gen-tournament-pop [pop : Population] [size : number]) : Population
  (cond
    [(= 0 size) empty]
    [else (cons (get-random pop)
                (gen-tournament-pop pop (sub1 size)))]))

; Gets a random member of a population.
(define (get-random [pop : Population]) : Tour
  (nth (add1 (random-range kPopulationSize)) pop))

;;; 4. Crossover
; Crosses over two tours by randomly selecting an element
; during iteration.
(define (crossover [i1 : Tour] [i2 : Tour]) : Tour
  (local ([define split (split-tour (floor (/ (length i1) 2)) i1)]
          [define left (fst split)]
          [define right (snd split)])
    (append left (reorder right i2))))

; Takes a list of cities and splits it into two over a split point
(define (split-tour [split : number] [tour : Tour])
  : (Tour * Tour)
  (cond
    [(= 0 split) (values empty tour)]
    [else (let ([recur (split-tour (sub1 split) (rest tour))])
            (values (cons (first tour) (fst recur))
                    (snd recur)))]))

; Takes a list of cities to find in a tour and returns them
; in the order in which they appear in the tour.
(define (reorder [cities : Tour] [tour : Tour])
  : Tour
  (cond
    [(empty? tour) empty]
    [(member (first tour) cities) (cons (first tour)
                                        (reorder cities (rest tour)))]
    [else (reorder cities (rest tour))]))

;;;; 5. Mutation
; Mutates a population.
(define (mutate-population [pop : Population]) : Population
  (map mutate-tour pop))

; Mutates an tour by randomly changing elements according to the
; mutation rate.
;  Mutation rate : kMutationRate
(define (mutate-tour [i : Tour]) : Tour
  (cond
    [(< (random) kMutationRate) (random-swap i)]
    [else i]))

; Randomly swaps two cities in a tour
(define (random-swap [tour : Tour]) : Tour
  (local ([define left (add1 (random-range (length tour)))]
          [define right (add1 (random-range (length tour)))])
    (swap-helper tour
               1
               left
               right
               (nth left tour)
               (nth right tour))))

; Helper for swapping
(define (swap-helper [tour : Tour]
                     [curr : number]
                     [index1 : number]
                     [index2 : number]
                     [val1 : City]
                     [val2 : City]) : Tour
  (cond
    [(empty? tour) empty]
    [(= curr index1)
     (cons val2 (swap-helper (rest tour) (add1 curr)
                             index1 index2 val1 val2))]
    [(= curr index2)
     (cons val1 (swap-helper (rest tour) (add1 curr)
                             index1 index2 val1 val2))]
    [else (cons (first tour) (swap-helper (rest tour) (add1 curr)
                                          index1 index2 val1 val2))]))

;;;;;;;;;;;;;;;;;;;;;;; ALGORITHM ;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Performs the genetic algorithm starting with a random population.
; Returns the best solution.
(define (genetic-algorithm) : void
  (local ([define world (generate-graph kGraphSize)]
          [define result (genetic-algorithm-helper (snd world)
                                                   (gen-pop (fst world)) 1)]
          [define generation (fst result)]
          [define solution (snd result)])
    (begin
      (displayln
       (string-append "Generations: "
                      (to-string generation)))
      (displayln
       (string-append "Solution: "
                      (to-string solution)))
      (displayln
       (string-append "Tour length: "
                      (to-string (calc-tour-length solution (snd world))))))))

; Recursive helper for genetic algorithm.
; Prints the generation and the current best solution according to flag.
; Performs iterations on population until the best solution is found.
; Returns the tour found and the number of generations produced.
(define (genetic-algorithm-helper [graph : Graph] [pop : Population]
                                  [generation : number])
  : (number * Tour)
  (local ([define best (get-best pop graph)]
          [define best-difference (calc-tour-length best graph)])
    (cond
      [(and (> kGenerationCap 0)
            (> generation kGenerationCap)) (values generation best)]
      [else
       (begin
         (if kPrintGenerations
             (displayln (list "Generation:" (to-string generation)
                              "  Best:"
                              (to-string (exact->inexact best-difference))))
             (void))
         (genetic-algorithm-helper graph (iteration pop graph)
                                   (add1 generation)))])))

; Iterates a population by performing tournament, crossover, and mutation.
(define (iteration [pop : Population] [graph : Graph]) : Population
  (cons (get-best pop graph)
          (mutate-population (iter-helper pop (sub1 kPopulationSize) graph))))

; Recursive helper for iteration.
(define (iter-helper [pop : Population] [n : number] [graph : Graph])
  : Population
  (cond
    [(= n 0) empty]
    [else (local ([define i1 (tournament pop graph)]
                  [define i2 (tournament pop graph)])
            (cons (crossover i1 i2) (iter-helper pop (sub1 n) graph)))]))

;;;;; TEST COVERAGE ;;;;;
(define backup-kPrintGenerations kPrintGenerations)

(set! kPrintGenerations true)
(genetic-algorithm)

(set! kPrintGenerations false)
(genetic-algorithm)

(set! kPrintGenerations backup-kPrintGenerations)