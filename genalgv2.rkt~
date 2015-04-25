#lang plai-typed
(print-only-errors true)
(require (typed-in racket
                   [build-list : (number (number -> number) -> (listof number))]
                   [displayln : ('a -> void)]
                   [exact->inexact : (number -> number)]))

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
; Secondary algorithm parameters include:
;  Individual size
;  Algorithm solution
;  Solution score
;  Acceptable fitness

; Secondary algorithm parameters are implementation
; dependent, as they are particular to the type of
; individuals in the population.

; To customize the algorithm, the secondary parameters
; must be configured, and the following functions must
; be re-implemented:
;  gen-individual
;  calc-fitness
;  crossover
;  mutate-individual
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

;;;;;;;;;; PARAMETERS ;;;;;;;;;;

;;;;; Types ;;;;;
(define-type-alias Individual (listof number))
(define-type-alias Population (listof Individual))

;;;;; Algorithm Parameters ;;;;;
; The size of an individual
(define kIndividualSize 64)

; The size of the population
(define kPopulationSize 50)

; The number of individuals selected for a tournament
(define kTournamentSize 5)

; The rate of random mutation
(define kMutationRate 0.015)

; Elitism flag
;  Elitism preserves the best individual of the population
;  for the next generation of the algorithm.
(define kElitism true)

; Generation cap
;  Caps the number of generations so that the best solution is
;  returned after a certain number of iterations.
; A generation cap of 0 means the algorithm is allowed to run
; until an acceptable solution is found.
(define kGenerationCap 0)

;;;;; Solution ;;;;;
(define kSolution (build-list kIndividualSize (lambda (x) 1)))

; Raw fitness score for solution.
(define kSolutionScore kIndividualSize)

; Minimum acceptable fitness for algorithm termination
(define kAcceptableFitness kSolutionScore)

;;;;;;;;;; ALGORITHM MACHINERY ;;;;;;;;;;

;;;; 1. INITIALIZATION ;;;;;
; Generates a population randomly.
;   Population size : kPopulationSize
(define (gen-pop) : Population
  (gen-pop-helper kPopulationSize))

; Generates an individual.
;  Indivual size : kIndividualSize
(define (gen-individual) : Individual
  (gen-individual-helper kIndividualSize))

; Recursive helper for population generation.
;  Generates a population of 
(define (gen-pop-helper [size : number]) : Population
  (cond
    [(= size 0) empty]
    [else (cons (gen-individual)
                (gen-pop-helper (sub1 size)))]))

; Recursive helper for individual generation.
(define (gen-individual-helper [size : number]) : Individual
  (cond
    [(= size 0) empty]
    [else (cons (random-range 2)
                (gen-individual-helper (sub1 size)))]))

;;;; 2. EVALUATION
; Evaluates the fitness of an individual.
;  Checks against : kSolution
(define (calc-fitness [indiv : Individual]) : number
  (fitness-helper indiv kSolution))

; Recursive helper for fitness evaluation.
(define (fitness-helper [indiv : Individual] [solution : Individual]) : number
  (cond
    [(empty? indiv) 0]
    [else (let ([recur (fitness-helper (rest indiv)
                                       (rest solution))])
            (if (= (first indiv) (first solution))
                (add1 recur)
                recur))]))

(test (calc-fitness (build-list kIndividualSize (lambda (x) 1))) 64)
(test (calc-fitness (build-list kIndividualSize (lambda (x) 0))) 0)

;;;; 3. SELECTION
; Gets the highest fitness member of a population.
(define (get-best [pop : Population]) : Individual
  (best-helper (rest pop) (first pop) (calc-fitness (first pop))))

; Recursive helper for best individual.
(define (best-helper [pop : Population]
                     [best : Individual]
                     [best-score : number]) : Individual
  (cond
    [(empty? pop) best]
    [else (local ([define curr (first pop)]
                  [define curr-score (calc-fitness curr)])
            (cond
              [(> curr-score best-score)
               (best-helper (rest pop) curr curr-score)]
              [else (best-helper (rest pop) best best-score)]))]))

; Performs a tournament and returns the winning individual.
;  Tournament size : kTournamentSize
(define (tournament [pop : Population]) : Individual
  (get-best (gen-tournament-pop pop kTournamentSize)))

; Generates a tournament population of the given size.
(define (gen-tournament-pop [pop : Population] [size : number]) : Population
  (cond
    [(= 0 size) empty]
    [else (cons (get-random pop)
                (gen-tournament-pop pop (sub1 size)))]))

; Gets a random member of a population.
(define (get-random [pop : Population]) : Individual
  (nth (add1 (random-range kPopulationSize)) pop))

;;; 4. Crossover
; Crosses over two individuals by randomly selecting an element
; during iteration.
(define (crossover [i1 : Individual] [i2 : Individual]) : Individual
  (cond
    [(empty? i1) empty]
    [else (cons
           (if (= 0 (random-range 2))
               (first i1)
               (first i2))
           (crossover (rest i1) (rest i2)))]))

;;;; 5. Mutation
; Mutates a population.
(define (mutate-population [pop : Population]) : Population
  (map mutate-individual pop))

; Mutates an individual by randomly changing elements according to the
; mutation rate.
;  Mutation rate : kMutationRate
(define (mutate-individual [i : Individual]) : Individual
  (cond
    [(empty? i) empty]
    [else (cons
           (if (< (random) kMutationRate)
               (random-range 2)
               (first i))
           (mutate-individual (rest i)))]))

;;;;;;;;;;;;;;;;;;;;;;; ALGORITHM ;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Performs the genetic algorithm starting with a random population.
; Returns the best solution.
(define (genetic-algorithm) : Individual
  (genetic-algorithm-helper (gen-pop) 1))

; Recursive helper for genetic algorithm.
; Prints the generation and the current best solution.
; Performs iterations on population until the best solution is found.
(define (genetic-algorithm-helper [pop : Population] [generation : number]) : Individual
  (local ([define best (get-best pop)]
          [define best-fitness (calc-fitness best)])
    (cond
      [(and (> kGenerationCap 0)
            (> generation kGenerationCap)) best]
      [else (begin
              (displayln (list "Generation:" (to-string generation)
                               "  Best:" (to-string (exact->inexact best-fitness))))
              (cond
                [(>= best-fitness kAcceptableFitness) best]
                [else (genetic-algorithm-helper (iteration pop) (add1 generation))]))])))

; Iterates a population by performing tournament, crossover, and mutation.
(define (iteration [pop : Population]) : Population
  (cond
    [kElitism (cons (get-best pop)
                   (mutate-population (iter-helper pop (sub1 kPopulationSize))))]
    [else (mutate-population (iter-helper pop kPopulationSize))]))

; Recursive helper for iteration.
(define (iter-helper [pop : Population] [n : number]) : Population
  (cond
    [(= n 0) empty]
    [else (local ([define i1 (tournament pop)]
                  [define i2 (tournament pop)])
            (cons (crossover i1 i2) (iter-helper pop (sub1 n))))]))

;;;;; TEST COVERAGE ;;;;;
(define backup-elitism kElitism)
(define backup-kGenerationCap kGenerationCap)

(set! kElitism false)
(set! kGenerationCap 0)
(genetic-algorithm)

(set! kElitism true)
(set! kGenerationCap 1)
(genetic-algorithm)

(set! kElitism backup-elitism)
(set! kGenerationCap backup-kGenerationCap)