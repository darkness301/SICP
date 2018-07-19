;#lang racket/base
;#lang racket
;(require racket/trace)
#lang planet neil/sicp

;;;;CODE FROM CHAPTER 1 OF STRUCTURE AND INTERPRETATION OF COMPUTER PROGRAMS

;;; Examples from the book are commented out with ;: so that they
;;;  are easy to find and so that they will be omitted if you evaluate a
;;;  chunk of the file (programs with intervening examples) in Scheme.

;;; BEWARE: Although the whole file can be loaded into Scheme,
;;;  don't expect the programs to work if you do so.  For example,
;;;  the redefinition of + in exercise 1.9 wreaks havoc with the
;;;  last version of square defined here.


;;;SECTION 1.1.1

;; interpreter examples

;: 486

;: (+ 137 349)
;: (- 1000 334)
;: (* 5 99)
;: (/ 10 5)
;: (+ 2.7 10)

;: (+ 21 35 12 7)
;: (* 25 4 12)

;: (+ (* 3 5) (- 10 6))

;: (+ (* 3 (+ (* 2 4) (+ 3 5))) (+ (- 10 7) 6))

;: (+ (* 3
;:       (+ (* 2 4)
;:          (+ 3 5)))
;:    (+ (- 10 7)
;:       6))


;;;SECTION 1.1.2

;: (define size 2)
;: size
;: (* 5 size)

;: (define pi 3.14159)
;: (define radius 10)
;: (* pi (* radius radius))
;: (define circumference (* 2 pi radius))
;: circumference


;;;SECTION 1.1.3
;
; (* (+ 2 (* 4 6))
;    (+ 3 5 7))


;;;SECTION 1.1.4

(define (square x) (* x x))

; (square 21)
; (square (+ 2 5))
; (square (square 3))

(define (sum-of-squares x y)
  (+ (square x) (square y)))

; (sum-of-squares 3 4)

;(define (f a)
;  (sum-of-squares (+ a 1) (* a 2)))

;: (f 5)


;;;SECTION 1.1.5

; Applicative oder

;: (f 5)
;: (sum-of-squares (+ 5 1) (* 5 2))
;: (+ (square 6) (square 10))
;: (+ (* 6 6) (* 10 10))
;: (+ 36 100)
;:  136

; Normal order

;: (f 5)
;: (sum-of-squares (+ 5 1) (* 5 2))
;: (+    (square (+ 5 1))      (square (* 5 2))  )
;: (+    (* (+ 5 1) (+ 5 1))   (* (* 5 2) (* 5 2)))
;: (+         (* 6 6)             (* 10 10))
;: (+           36                   100)
;:                     136


;;;SECTION 1.1.6

;(define (abs x)
;  (cond ((> x 0) x)
;        ((= x 0) 0)
;        ((< x 0) (- x))))

(define (abs x)
  (cond ((< x 0) (- x))
        (else x)))

;(define (abs x)
;  (if (< x 0)
;      (- x)
;      x))

;: (and (> x 5) (< x 10))

(define (>= x y)
  (or (> x y) (= x y)))

;(define (>= x y)
;  (not (< x y)))


;;EXERCISE 1.1

; 10
;
; (+ 5 3 4)
;
; (- 9 1)
;
; (/ 6 2)
;
; (+ (* 2 4) (- 4 6))

; (define a 3)
;
; (define b (+ a 1))
;
; (+ a b (* a b))
;
; (= a b)
;
; (if (and (> b a) (< b (* a b)))
;     b
;     a)
;
; (cond ((= a 4) 6)
;       ((= b 4) (+ 6 7 a))
;       (else 25))
;
; (+ 2 (if (> b a) b a))
;
; (* (cond ((> a b) a)
; 	 ((< a b) b)
; 	 (else -1))
;    (+ a 1))

;;EXERCISE 1.2

;(/ (+ 5 4 (- 2 (- 3 (+ 6 (/ 4 5)))))
;      (* 3 (- 6 2) (- 2 7)))

;;EXERCISE 1.3

;(define (bigger-sum-of-squares x y z)
;  (cond ((and (>= x y) (>= y z)) (sum-of-squares x y))
;        ((and (>= x y) (>= z y)) (sum-of-squares x y))
;        ((and (>= y x) (>= x z)) (sum-of-squares y x))
;        ((and (>= y x) (>= z x)) (sum-of-squares x y))))

(define (bigger x y)
  (if (> x y)
      x
      y))

(define (smaller x y)
  (if (> x y)
      y
      x))

(define (bigger-sum-of-squares x y z)
  (sum-of-squares (bigger x y)                    
                  (bigger (smaller x y) z))) 

;;EXERCISE 1.4

;(define (a-plus-abs-b a b)
;  ((if (> b 0) + -) a b))

;;EXERCISE 1.5

;(define (p) (p))
;
;(define (test x y)
;  (if (= x 0)
;      0
;      y))

;: (test 0 (p))


;;;SECTION 1.1.7

;(define (sqrt-iter guess x)
;  (if (good-enough? guess x)
;      guess
;      (sqrt-iter (improve guess x)
;                 x)))

;(define (improve guess x)
;  (average guess (/ x guess)))

(define (average x y)
  (/ (+ x y) 2))

;(define (good-enough? guess x)
;  (< (abs (- (square guess) x)) 0.001))

;(define (sqrt x)
;  (sqrt-iter 1.0 x))


; (sqrt 9)
; (sqrt (+ 100 37))
; (sqrt (+ (sqrt 2) (sqrt 3)))
; (square (sqrt 1000))
;
;
;;;EXERCISE 1.6

(define (new-if predicate then-clause else-clause)
  (cond (predicate then-clause)
        (else else-clause)))
;
; (new-if (= 2 3) 0 5)
;
; (new-if (= 1 1) 0 5)

;(define (sqrt-iter guess x)
;  (new-if (good-enough? guess x)
;          guess
;          (sqrt-iter (improve guess x)
;                     x)))

;(trace sqrt-iter)

;;; Test "if" and "new-if":

; (if #t (display "good") (display "bad"))
; (new-if #t (display "good") (display "bad"))

;;; EXERCISE 1.7

;(define (good-enough? old-guess new-guess)
;    (> 0.01
;       (/ (abs (- new-guess old-guess))
;          old-guess)))
;
;(define (sqrt-iter guess x)
;    (if (good-enough? guess (improve guess x))  
;        (improve guess x)
;        (sqrt-iter (improve guess x)
;                   x)))

;(define (sqrt-iter old-guess x)
;    (let ((new-guess (improve old-guess x)))
;        (if (good-enough? old-guess new-guess)
;            new-guess
;            (sqrt-iter new-guess x))))

;;; EXERCISE 1.8

;(define (cube x)
;  (* x x x))
;
;(define (cube-root x)
;  (cube-root-iter 1.0 x))
;
;(define (cube-root-iter guess x)           
;  (if (good-enough? guess x)              
;      guess
;      (cube-root-iter (improve guess x)
;                      x)))
;
;(define (good-enough? guess x)             
;  (< (abs (- (cube guess) x))
;     0.001))
;
;(define (improve guess x)                   
;  (/ (+ (/ x (square guess)) (* 2 guess))
;     3))

;;; SECTION 1.1.8

;(define (square x) (* x x))

;(define (square x) 
;  (exp (double (log x))))

;(define (double x) (+ x x))
;
;
;;; As in 1.1.7

;(define (sqrt x)
;  (sqrt-iter 1.0 x))
;
;(define (sqrt-iter guess x)
;  (if (good-enough? guess x)
;      guess
;      (sqrt-iter (improve guess x) x)))
;
;(define (good-enough? guess x)
;  (< (abs (- (square guess) x)) 0.001))
;
;(define (improve guess x)
;  (average guess (/ x guess)))


;;; Block-structured

;(define (sqrt x)
;  (define (good-enough? guess x)
;    (< (abs (- (square guess) x)) 0.001))
;  (define (improve guess x)
;    (average guess (/ x guess)))
;  (define (sqrt-iter guess x)
;    (if (good-enough? guess x)
;        guess
;        (sqrt-iter (improve guess x) x)))
;  (sqrt-iter 1.0 x))

;;; Taking advantage of lexical scoping

;(define (sqrt x)
;  (define (good-enough? guess)
;    (< (abs (- (square guess) x)) 0.001))
;  (define (improve guess)
;    (average guess (/ x guess)))
;  (define (sqrt-iter guess)
;    (if (good-enough? guess)
;        guess
;        (sqrt-iter (improve guess))))
;  (sqrt-iter 1.0))

;;;;SECTION 1.2.1
;
;;; Recursive

;(define (factorial n)
;  (if (= n 1)
;      1
;      (* n (factorial (- n 1)))))

;;; Iterative

;(define (factorial n)
;  (fact-iter 1 1 n))
;
;(define (fact-iter product counter max-count)
;  (if (> counter max-count)
;      product
;      (fact-iter (* counter product)
;                 (+ counter 1)
;                 max-count)))


;;; Iterative, block-structured (from footnote)

;(define (factorial n)
;  (define (iter product counter)
;    (if (> counter n)
;        product
;        (iter (* counter product)
;              (+ counter 1))))
;  (iter 1 1))

;;;EXERCISE 1.9

;(define (inc x)
;  (+ x 1))
;
;(define (dec x)
;  (- x 1))

;(define (＋ a b)
;  (if (= a 0)
;      b
;      (inc (＋ (dec a) b))))

;(define (＋ a b)
;  (if (= a 0)
;      b
;      (＋ (dec a) (inc b))))

;;;EXERCISE 1.10

;(define (A x y)
;  (cond ((= y 0) 0)
;        ((= x 0) (* 2 y))
;        ((= y 1) 2)
;        (else (A (- x 1)
;                 (A x (- y 1))))))
;
;(A 1 10)
;(A 2 4)
;(A 3 3)
;
;;2n
;(define (f n) (A 0 n))

;;2 exp n
;(define (g n) (A 1 n))

;;2 exp 2 exp 2....
;(define (h n) (A 2 n))
;
;(define (k n) (* 5 n n))

;;;;SECTION 1.2.2
;
;;; Recursive

;(define (fib n)
;  (cond ((= n 0) 0)
;        ((= n 1) 1)
;        (else (+ (fib (- n 1))
;                 (fib (- n 2))))))

;;; Iterative

;(define (fib n)
;  (fib-iter 1 0 n))
;
;(define (fib-iter a b count)
;  (if (= count 0)
;      b
;      (fib-iter (+ a b) a (- count 1))))


;;; Counting change

(define (count-change amount)
  (cc amount 5))

(define (cc amount kinds-of-coins)
  (cond ((= amount 0) 1)
        ((or (< amount 0) (= kinds-of-coins 0)) 0)
        (else (+ (cc amount
                     (- kinds-of-coins 1))
                 (cc (- amount
                        (first-denomination kinds-of-coins))
                     kinds-of-coins)))))

(define (first-denomination kinds-of-coins)
  (cond ((= kinds-of-coins 1) 1)
        ((= kinds-of-coins 2) 5)
        ((= kinds-of-coins 3) 10)
        ((= kinds-of-coins 4) 25)
        ((= kinds-of-coins 5) 50)))

;(count-change 100)

;;EXERCISE 1.11

;rec-version

;(define (f n)
;    (if (< n 3)
;        n
;        (+ (f (- n 1))
;           (* 2 (f (- n 2)))
;           (* 3 (f (- n 3))))))

;iter-version

;(define (f n)
;  (f-iter 2 1 0 0 n))
;
;(define (f-iter a b c i n)
;  (if (= i n)
;      c
;      (f-iter (+ a (* 2 b) (* 3 c))  
;              a                      
;              b                       
;              (+ i 1)
;              n)))

;;EXERCISE 1.12

; Recursive Version

;(define (pascal row col)
;    (cond ((> col row)
;            (error "invalid col value"))
;          ((or (= col 0) (= row col))
;            1)
;          (else (+ (pascal (- row 1) (- col 1))
;                   (pascal (- row 1) col)))))

; Iterative Version

(define (pascal row col)
  (/ (factorial row)
     (* (factorial col)
        (factorial (- row col)))))

;;;;SECTION 1.2.3

;;;EXERCISE 1.14

;(count-change 11)
;(cc 11 5)
;(+ (cc 11 4) (cc (- 11 50) 5))
;(+ (+ (cc 11 3) (cc (- 11 25) 4)) (cc -39 5))
;(+ (+ (+ (cc 11 2) (cc (- 11 10) 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (cc 11 1) (cc (- 11 5) 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ (cc 11 0) (cc (- 11 1) 1)) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ (cc 11 0) (cc 10 1)) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (cc 10 1)) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ (cc 10 0) (cc (- 10 1) 1))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (cc 9 1))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ (cc 9 0) (cc (- 9 1) 1)))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ 0 (cc 8 1)))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ 0 (+ (cc 8 0) (cc (- 8 1) 1))))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ 0 (+ 0 (cc 7 1))))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ 0 (+ 0 (+ (cc 7 0) (cc (- 7 1) 1)))))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ 0 (+ 0 (+ 0 (cc 6 1)))))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ 0 (+ 0 (+ 0 (+ (cc 6 0) (cc (- 6 1) 1))))))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ 0 (+ 0 (+ 0 (+ 0 (cc 5 1))))))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;(+ (+ (+ (+ (+ 0 (+ 0 (+ 0 (+ 0 (+ 0 (+ 0 (+ (cc 5 0) (cc (- 5 1) 1)))))))) (cc 6 2)) (cc 1 3)) (cc -14 4)) (cc -39 5))
;; etc...

; so (cc n 1) generates 2*n steps, with n space taken up
;    O(n) in both time and space
; (cc n 2) generates a (cc n 1) process and a (cc (-n 1) 2) process, which means
;    looks recursive, right? so for every increase of the number of coins, we're
;    going to have a new branch of O(n) to worry about, so it's O(n^5) in time and space

;;;EXERCISE 1.16

(define (cube x) (* x x x))

(define (p x) (- (* 3 x) (* 4 (cube x))))

(define (sine angle)
  (if (not (> (abs angle) 0.1))
      angle
      (p (sine (/ angle 3.0)))))

;(sine 12.15)
;(p (sine (/ 12.15 3.0)))
;(p (sine 4.05))
;(p (p (sine (/ 4.05 3.0))))
;(p (p (sine 1.3499999999999999)))
;(p (p (p (sine (/ 1.3499999999999999 3.0)))))
;(p (p (p (sine 0.44999999999999996))))
;(p (p (p (p (sine (/ 0.44999999999999996 3.0))))))
;(p (p (p (p (sine 0.15)))))
;(p (p (p (p (p (sine (/ 0.15 3.0)))))))
;(p (p (p (p (p (sine 4.9999999999999996e-2))))))
;(p (p (p (p (p 4.9999999999999996e-2)))))

;; part a) p gets applied 5 times

;; part b) O(log n) in both time and space (dividing by 3 each time)

;;;;SECTION 1.2.4

;;; Linear recursion

;(define (expt b n)
;  (if (= n 0)
;      1
;      (* b (expt b (- n 1)))))

;;; Linear iteration

;(define (expt b n)
;  (expt-iter b n 1))
;
;(define (expt-iter b counter product)
;  (if (= counter 0)
;      product
;      (expt-iter b
;                (- counter 1)
;                (* b product)))) 

;;; Logarithmic iteration

;(define (fast-expt b n)
;  (cond ((= n 0) 1)
;        ((even? n) (square (fast-expt b (/ n 2))))
;        (else (* b (fast-expt b (- n 1))))))
;
;(define (even? n)
;  (= (remainder n 2) 0))

;EXERCISE 1.16

(define (fast-expt b n)
  (expt-iter b n 1))

(define (expt-iter b n a)
  (cond ((= n 0)
         a)
        ((even? n)
         (expt-iter (square b)
                    (/ n 2)
                    a))
        ((odd? n)
         (expt-iter b
                    (- n 1)
                    (* b a)))))

;;;EXERCISE 1.17

;(define (* a b)
;  (if (= b 0)
;      0
;      (+ a (* a (- b 1)))))

;(define (fast-multiply a b)
;  (define (double x) (* x 2))
;  (define (halve x) (/ x 2))
;
;  (cond ((or (= 0 a) (= 0 b)) 0)
;        ((even? b) (fast-multiply (double a) (halve b)))
;        (else (+ a (fast-multiply a (- b 1))))))


;;;EXERCISE 1.18

(define (fast-multiply a b)
  (define (double x) (* x 2))
  (define (halve x) (/ x 2))
  
  (define (iter a b acc)
    (cond ((< b 1) acc)
          ((even? b) (iter (double a) (halve b) acc))
          (else (iter a (- b 1) (+ a acc)))))
  
  (iter a b 0))

;;;EXERCISE 1.19

(define (fib n)
  (fib-iter 1 0 0 1 n))

(define (fib-iter a b p q count)
  (cond ((= count 0) b)
        ((even? count)
         (fib-iter a
                   b
                   (+ (square p) (square q))     
                   (+ (* 2 p q) (square q))      
                   (/ count 2)))
        (else (fib-iter (+ (* b q) (* a q) (* a p))
                        (+ (* b p) (* a q))
                        p
                        q
                        (- count 1)))))


;;;;SECTION 1.2.5

(define (gcd a b)
  (if (= b 0)
      a
      (gcd b (remainder a b))))

;(trace gcd)

;;;EXERCISE 1.20

; Normal order:

;(gcd 206 40)
;(if (= 40 0)
;    206
;    (gcd 40 (remainder 206 40)))
;(gcd 40 (remainder 206 40))
;(if (= (remainder 206 40) 0)
;    40
;     (gcd (remainder 206 40) (remainder 40 (remainder 206 40))))
; ; remainder evaluations = 1
;(if (= 6 0)
;    40
;    (gcd (remainder 206 40) (remainder 40 (remainder 206 40))))
;(gcd (remainder 206 40) (remainder 40 (remainder 206 40)))
;(if (= (remainder 40 (remainder 206 40)) 0)
;    (remainder 206 40)
;    (gcd (remainder 40 (remainder 206 40)) (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))
; ; remainder evaluations = 2
;(if (= (remainder 40 6) 0)
;    (remainder 206 40)
;    (gcd (remainder 40 (remainder 206 40)) (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))
; ; remainder evaluations = 3
;(if (= 4 0)
;    (remainder 206 40)
;    (gcd (remainder 40 (remainder 206 40)) (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))
;(gcd (remainder 40 (remainder 206 40))
;     (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))
;(if (= (remainder (remainder 206 40) (remainder 40 (remainder 206 40))) 0)
;    (remainder 40 (remainder 206 40))
;    (gcd (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;         (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))))
; ; remainder evaluations = 5
;(if (= (remainder 6 (remainder 40 6)) 0)
;    (remainder 40 (remainder 206 40))
;    (gcd (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;         (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))))
; ; remainder evaluations = 6
;(if (= (remainder 6 4) 0)
;    (remainder 40 (remainder 206 40))
;    (gcd (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;         (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))))
; ; remainder evaluations = 7
;(if (= 2 0)
;    (remainder 40 (remainder 206 40))
;    (gcd (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;         (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))))
;
;(gcd (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;     (remainder (remainder 40 (remainder 206 40))
;                (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))
;
;(if (= (remainder (remainder 40 (remainder 206 40))
;                  (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))
;       0)
;    (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;    (gcd (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))
;         (remainder (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;                    (remainder (remainder 40 (remainder 206 40))
;                               (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))))
;
; ; remainder evaluated: 10
;(if (= (remainder (remainder 40 6)
;                  (remainder 6 (remainder 40 6)))
;       0)
;    (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;    (gcd (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))
;         (remainder (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;                    (remainder (remainder 40 (remainder 206 40))
;                               (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))))
; ; remainder evaluated: 12
;(if (= (remainder 4 (remainder 6 4)) 0)
;    (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;    (gcd (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))
;         (remainder (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;                    (remainder (remainder 40 (remainder 206 40))
;                               (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))))
; ; remainder evaluated: 13
;(if (= (remainder 4 2) 0)
;    (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;    (gcd (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))
;         (remainder (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;                    (remainder (remainder 40 (remainder 206 40))
;                               (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))))
; ; remainder evaluated: 14
;(if (= 0 0)
;    (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;    (gcd (remainder (remainder 40 (remainder 206 40))
;                    (remainder (remainder 206 40) (remainder 40 (remainder 206 40))))
;         (remainder (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
;                    (remainder (remainder 40 (remainder 206 40))
;                               (remainder (remainder 206 40) (remainder 40 (remainder 206 40)))))))
;(remainder (remainder 206 40) (remainder 40 (remainder 206 40)))
; ; remainder evaluated: 16
;(remainder 6 (remainder 40 6))
; ; remainder evaluated: 17
;(remainder 6 4)
; ; remainder evaluated: 18
;2
;
;;;; WHEW! Normal order took 18 applications of (remainder a b)
;
;
;; Applicative:
;
;(gcd 206 40)
;(if (= 40 0)
;    206
;    (gcd 40 (remainder 206 40)))
;(if #f 206 (gcd 40 (remainder 206 40)))
;; 1 remainder
;(gcd 40 (remainder 206 40))
;(gcd 40 6)
;(if (= 6 0)
;    40
;    (gcd 6 (remainder 40 6)))
;(if #f 40 (gcd 6 (remainder 40 6)))
;(gcd 6 (remainder 40 6))
;; 2 remainder
;(gcd 6 4)
;(if (= 4 0)
;    6
;    (gcd 4 (remainder 6 4)))
;
;(if #f 6 (gcd 4 (remainder 6 4)))
;(gcd 4 (remainder 6 4))
;; 3 remainder
;(gcd 4 2)
;(if (= 2 0)
;    4
;    (gcd 2 (remainder 4 2)))
;(if #f 4 (gcd 2 (remainder 4 2)))
;(gcd 2 (remainder 4 2))
;; 4 remainder
;(gcd 2 0)
;(if (= 0 0)
;    2
;    (gcd 0 (remainder 2 0)))
;(if #t 2 (gcd 0 (remainder 2 0)))
;2
;
;;;; That was easier: applicative order takes 4 applications of (remainder a b)

;;;;SECTION 1.2.6
;
;;; prime?

(define (smallest-divisor n)
  (find-divisor n 2))

; base the fact that if d is a disvisor of n, then so is n/d.But d and n/d cannot both be greater than sqrt of n.

;(define (find-divisor n test-divisor)
;  (cond ((> (square test-divisor) n) n)
;        ((divides? test-divisor n) test-divisor)
;        (else (find-divisor n (+ test-divisor 1)))))

(define (divides? a b)
  (= (remainder b a) 0))

;(define (prime? n)
;  (= n (smallest-divisor n)))

;; fast-prime?

;(define (expmod base exp m)
;  (cond ((= exp 0) 1)
;        ((even? exp)
;         (remainder (square (expmod base (/ exp 2) m))
;                    m))
;        (else
;         (remainder (* base (expmod base (- exp 1) m))
;                    m))))        

(define (fermat-test n)
  (define (try-it a)
    (= (expmod a n n) a))
  (try-it (+ 1 (random (- n 1)))))

(define (fast-prime? n times)
  (cond ((= times 0) true)
        ((fermat-test n) (fast-prime? n (- times 1)))
        (else false)))

;;;EXERCISE 1.21

;(smallest-divisor 199)
;(smallest-divisor 1999)
;(smallest-divisor 19999)

;;;EXERCISE 1.22

;(define (runtime) 
;  (current-milliseconds))

;(define (even? n)
;  (= 0 (remainder n 2)))

(define (timed-prime-test n)
  (start-prime-test n (runtime)))

(define (start-prime-test n start-time)
  (if (prime? n)
      (report-prime n (- (runtime) start-time))
      #f))

(define (report-prime n elapsed-time)
  (display "Prime is ")
  (display n)
  (display " and time elasped ")
  (display elapsed-time)
  (display "ms\n"))

;search counter 
(define (search-for-primes a n)
  (search-for-primes-helper (next-odd a) 0 n))

(define (search-for-primes-helper a found n)
  (if (= found n)
      (display "Done.")
      (search-for-primes-helper (+ a 2)
                                (if (timed-prime-test a)
                                    (+ found 1)
                                    found)
                                n)))
(define (next-odd n)
  (if (even? n)
      (+ n 1)
      n))


;;;EXERCISE 1.23

(define (next n)
  (if (= n 2)
      3
      (+ n 2)))

(define (find-divisor n test-divisor)
  (cond ((> (square test-divisor) n) 
         n)
        ((divides? test-divisor n) 
         test-divisor)
        (else 
         (find-divisor n (next test-divisor)))))

;;;EXERCISE 1.24

(define (prime? n)
  (fast-prime? n 10))

;EXERCISE 1.25

;(define (expmod base exp m)
;  (remainder (fast-expt base exp) m))

;;EXERCISE 1.26

;(define (expmod base exp m)
;  (cond ((= exp 0) 1)
;        ((even? exp)
;         (remainder (* (expmod base (/ exp 2) m)
;                       (expmod base (/ exp 2) m))
;                    m))
;        (else
;         (remainder (* base (expmod base (- exp 1) m))
;                    m))))


;;EXERCISE 1.27

(define (carmichael-test n)
  (test-iter 1 n))

;(define (test-iter a n)
;    (cond ((= a n)
;            #t)
;          ((congruent? a n)
;            (test-iter (+ a 1) n))
;          (else
;            #f)))

(define (congruent? a n)           
  (= (expmod a n n) a))

;test 
;(carmichael-test 561)
;(carmichael-test 1105)
;(carmichael-test 1729)
;(carmichael-test 2465)
;(carmichael-test 2821)
;(carmichael-test 6601)

;;EXERCISE 1.28

(define (expmod base exp m)
  (cond ((= exp 0)
         1)
        ((nontrivial-square-root? base m)                 
         0)                                              
        ((even? exp)
         (remainder (square (expmod base (/ exp 2) m))
                    m))
        (else
         (remainder (* base (expmod base (- exp 1) m))
                    m))))

(define (nontrivial-square-root? a n)
  (and (not (= a 1))
       (not (= a (- n 1)))
       (= 1 (remainder (square a) n))))

(define (non-zero-random n)
  (let ((r (random n)))
    (if (not (= r 0))
        r
        (non-zero-random n))))

(define (test-iter n times)
  (cond ((= times 0)
         #t)
        ((= (expmod (non-zero-random n) (- n 1) n)
            1)
         (test-iter n (- times 1)))
        (else
         #f)))

; Miller-Rabin-test display: undefined; miller-rabin-test work well, why?

(define (miller-rabin-test n)
  (let ((times (ceiling (/ n 2))))
    (test-iter n times)))

;test
;
;(Miller-Rabin-test 561)
;(Miller-Rabin-test 1105)
;(Miller-Rabin-test 1729)
;(Miller-Rabin-test 2465)
;(Miller-Rabin-test 2821)
;(Miller-Rabin-test 6601)


;;;;SECTION 1.3

;(define (cube x) (* x x x))

;;;;SECTION 1.3.1

;(define (sum-integers a b)
;  (if (> a b)
;      0
;      (+ a (sum-integers (+ a 1) b))))

;(define (sum-cubes a b)
;  (if (> a b)
;      0
;      (+ (cube a) (sum-cubes (+ a 1) b))))

;(define (pi-sum a b)
;  (if (> a b)
;      0
;      (+ (/ 1.0 (* a (+ a 2))) (pi-sum (+ a 4) b))))

;(define (sum term a next b)
;  (if (> a b)
;      0
;      (+ (term a)
;         (sum term (next a) next b))))

;;; Using sum

(define (inc n) (+ n 1))

(define (sum-cubes a b)
  (sum cube a inc b))

;(sum-cubes 1 10)


(define (identity x) x)

(define (sum-integers a b)
  (sum identity a inc b))

;(sum-integers 1 10)


(define (pi-sum a b)
  (define (pi-term x)
    (/ 1.0 (* x (+ x 2))))
  (define (pi-next x)
    (+ x 4))
  (sum pi-term a pi-next b))

;(* 8 (pi-sum 1 1000))

; Compute integral of a function

(define (integral f a b dx)
  (define (add-dx x) (+ x dx))
  (* (sum f (+ a (/ dx 2)) add-dx b)
     dx))

;(integral cube 0 1 0.01)
;(integral cube 0 1 0.001)

;;;EXERCISE 1.29

(define (simpson f a b n)
  (define h (/ (- b a) n))
  (define (y k)
    (f (+ a (* k h))))
  (define (factor k)
    (cond ((or (= k 0) (= k n))
           1)
          ((odd? k)
           4)
          (else
           2)))
  (define (term k)
    (* (factor k)
       (y k)))
  (define (next k)
    (+ k 1))
  (if (not (even? n))
      (error "n can't be odd")
      (* (/ h 3)
         (sum term (exact->inexact 0) next n))))

;(define (integral f a b)
;(define n 1000)
;(define h (/ (- b a) n)) 
;(define (add-2dx x) (+ x h h))
;(define h/3 (/ h 3)) 
;
;(* h/3 (+ (f a) 
;(* 2 (sum f (+ a h h) add-2dx (- b h h)))
;(* 4 (sum f (+ a h) add-2dx (- b h)))
;(f b))))

;;;EXERCISE 1.30

;(define (sum term a next b)
;  (define (iter a result)
;    (if (> a b)
;        result
;        (iter (next a)
;              (+ (term a) result))))
;  (iter a 0))

;;;EXERCISE 1.31

;(define (product term a next b)
;    (if (> a b)
;        1
;        (* (term a)
;           (product term (next a) next b))))

;(define (product term a next b)
;  (define (iter a result)
;    (if (> a b)
;        result
;        (iter (next a)
;              (* (term a) result))))
;  (iter a 1))

(define (factorial n)
  (product (lambda (x) x)
           1
           (lambda (i) (+ i 1))
           n))

(define (numer-term i)
  (cond ((= i 1)
         2)
        ((even? i)
         (+ i 2))
        (else
         (+ i 1))))

(define (denom-term i)
  (if (odd? i)
      (+ i 2)
      (+ i 1)))

(define (pi n)
  (* 4
     (exact->inexact
      (/ (product numer-term
                  1
                  (lambda (i) (+ i 1))
                  n)
         (product denom-term 
                  1
                  (lambda (i) (+ i 1))
                  n)))))

;;;EXERCISE 1.32

;(define (accumulate combiner null-value term a next b)
;  (if (> a b)
;      null-value
;      (combiner (term a)
;                (accumulate combiner
;                            null-value
;                            term
;                            (next a)
;                            next
;                            b))))

(define (sum term a next b)
  (accumulate + 
              0 
              term 
              a 
              next 
              b))

(define (product term a next b)
  (accumulate *
              1 
              term
              a
              next
              b))

(define (accumulate combiner null-value term a next b)
  (define (iter a result)
    (if (> a b)
        result
        (iter (next a)
              (combiner result (term a)))))
  (iter a null-value))

;;;EXERCISE 1.33

;(define (filtered-accumulate combine null-value term a next b valid?)
;  (if (> a b)
;      null-value
;      (let ((rest-terms (filtered-accumulate combine
;                                             null-value
;                                             term
;                                             (next a)
;                                             next
;                                             b
;                                             valid?)))
;        (if (valid? a)
;            (combine (term a) rest-terms)
;            rest-terms))))

(define (primes-sum a b)
  (filtered-accumulate + 
                       0
                       (lambda (x) x)
                       a
                       (lambda (i) (+ i 1))
                       b
                       prime?))

(define (coprime? i n)
  (and (< i n)
       (= 1 (gcd i n))))

(define (product-of-coprimes n)
  (filtered-accumulate *
                       1
                       (lambda (x) x)
                       1
                       (lambda (i) (+ i 1))
                       n
                       (lambda (x) (coprime? x n))))

(define (filtered-accumulate combine null-value term a next b valid?)
  (define (iter i result)
    (cond ((> i b)
           result)
          ((valid? i)
           (iter (next i) (combine (term i) result)))
          (else 
           (iter (next i) result))))
  (iter a null-value))

;;;;SECTION 1.3.2

;(define (pi-sum a b)
;  (sum (lambda (x) (/ 1.0 (* x (+ x 2))))
;       a
;       (lambda (x) (+ x 4))
;       b))
;
;(define (integral f a b dx)
;  (* (sum f
;          (+ a (/ dx 2.0))
;          (lambda (x) (+ x dx))
;          b)
;     dx))

;(define (plus4 x) (+ x 4))
;
;(define plus4 (lambda (x) (+ x 4)))

;((lambda (x y z) (+ x y (square z))) 1 2 3)
;
;
;;; Using let
;
;(define (f x y)
;  (define (f-helper a b)
;    (+ (* x (square a))
;       (* y b)
;       (* a b)))
;  (f-helper (+ 1 (* x y)) 
;            (- 1 y)))

;(define (f x y)
;  ((lambda (a b)
;     (+ (* x (square a))
;        (* y b)
;        (* a b)))
;   (+ 1 (* x y))
;   (- 1 y)))

;(define (f x y)
;  (let ((a (+ 1 (* x y)))
;        (b (- 1 y)))
;    (+ (* x (square a))
;       (* y b)
;       (* a b))))
;
;(+ (let ((x 3))
;     (+ x (* x 10)))
;   5)

;(let ((x 3)
;      (y (+ x 2)))
;  (* x y))

;(define (f x y)
;  (define a (+ 1 (* x y)))
;  (define b (- 1 y))
;  (+ (* x (square a))
;     (* y b)
;     (* a b)))
;
;
;;;EXERCISE 1.34

(define (f g)
  (g 2))

;(require racket/trace)
;(trace f)

;(f square)
;(f (lambda (z) (* z (+ z 1))))

;(f f)
;
;(f (lambda (g) (g 2)))
;
;((lambda (g) (g 2)) (lambda (g) (g 2)))
;
;((lambda (g) (g 2)) 2)
;
;(2 2)

;;;SECTION 1.3.3
;
;;; Half-interval method

(define (search f neg-point pos-point)
  (let ((midpoint (average neg-point pos-point)))
    (if (close-enough? neg-point pos-point)
        midpoint
        (let ((test-value (f midpoint)))
          (cond ((positive? test-value)
                 (search f neg-point midpoint))
                ((negative? test-value)
                 (search f midpoint pos-point))
                (else midpoint))))))

(define (close-enough? x y)
  (< (abs (- x y)) 0.001))

(define (half-interval-method f a b)
  (let ((a-value (f a))
        (b-value (f b)))
    (cond ((and (negative? a-value) (positive? b-value))
           (search f a b))
          ((and (negative? b-value) (positive? a-value))
           (search f b a))
          (else
           (error "Values are not of opposite sign" a b)))))


;(half-interval-method sin 2.0 4.0)
;
;(half-interval-method (lambda (x) (- (* x x x) (* 2 x) 3))
;                      1.0
;                      2.0)


;; Fixed points

(define tolerance 0.00001)

(define (fixed-point f first-guess)
  (define (close-enough? v1 v2)
    (< (abs (- v1 v2)) tolerance))
  (define (try guess)
    (let ((next (f guess)))
      (if (close-enough? guess next)
          next
          (try next))))
  (try first-guess))

;(fixed-point cos 1.0)
;
;(fixed-point (lambda (y) (+ (sin y) (cos y)))
;             1.0)


;(define (sqrt x)
;  (fixed-point (lambda (y) (/ x y))
;               1.0))
;
;(define (sqrt x)
;  (fixed-point (lambda (y) (average y (/ x y)))
;               1.0))


;;EXERCISE 1.35

;(define golden-ratio
;    (fixed-point (lambda (x) 
;                     (+ 1 (/ 1 x)))
;                 1.0))

;;EXERCISE 1.36
;
;(define tolerance 0.000001)
;
;(define (fixed-point f first-guess)
;  (define (close-enough? v1 v2)
;    (< (abs (- v1 v2)) tolerance))
;  (define (try guess step)
;    (display-info guess step)                       
;    (let ((next (f guess)))
;      (if (close-enough? next guess)
;          (begin                                  
;            (display-info next (+ 1 step))      
;            next)
;          (try next (+ 1 step)))))
;  (try first-guess 1))
;
;(define (display-info guess step)
;  (display "Step: ")
;  (display step)
;  (display " ")
;  
;  (display "Guess: ")
;  (display guess)
;  (newline))
;
;(define formula 
;  (lambda (x)
;    (/ (log 1000) 
;       (log x))))

;Test

;(fixed-point formula 2.0)
;(fixed-point (average-damp formula) 2.0)


;;;EXERCISE 1.37

;(cont-frac (lambda (i) 1.0)
;           (lambda (i) 1.0)
;           k)


;(define (cont-frac N D k)
;  (define (cf i)
;    (if (= k i)
;        (/ (N k) (D k))
;        (/ (N i)
;           (+ (D i) (cf (+ i 1))))))
;  (cf 1))

(define (cont-frac N D k)
  (define (iter i result)
    (if (= i 0)
        result
        (iter (- i 1)
              (/ (N i)
                 (+ (D i) result)))))
  (iter (- k 1)
        (/ (N k) (D k))))


(define (golden-ratio k)
  (+ 1
     (cont-frac (lambda (i) 1.0)
                (lambda (i) 1.0)
                k)))


;;;EXERCISE 1.38

(define (e k)
  (define (N i) 1)
  (define (D i)
    (if (= 0 (remainder (+ i 1) 3))
        (* 2 (/ (+ i 1) 3))
        1))
  (+ 2.0 
     (cont-frac N D k)))


;;;EXERCISE 1.39

(define (tan-cf x k)
  (define (N i)
    (if (= i 1)
        x
        (- (square x))))
  (define (D i)
    (- (* i 2) 1))
  (exact->inexact (cont-frac N D k)))


;;;;SECTION 1.3.4

(define (average-damp f)
  (lambda (x) (average x (f x))))

;((average-damp square) 10)

;(define (sqrt x)
;  (fixed-point (average-damp (lambda (y) (/ x y)))
;               1.0))
;
;(define (cube-root x)
;  (fixed-point (average-damp (lambda (y) (/ x (square y))))
;               1.0))


;;; Newton's method

(define (deriv g)
  (lambda (x)
    (/ (- (g (+ x dx)) (g x))
       dx)))

(define dx 0.00001)

;(define (cube x) (* x x x))
;((deriv cube) 5)

(define (newton-transform g)
  (lambda (x)
    (- x (/ (g x) ((deriv g) x)))))

(define (newtons-method g guess)
  (fixed-point (newton-transform g) guess))

;(define (sqrt x)
;  (newtons-method (lambda (y) (- (square y) x))
;                  1.0))

;;; Fixed point of transformed function

(define (fixed-point-of-transform g transform guess)
  (fixed-point (transform g) guess))

(define (sqrt x)
  (fixed-point-of-transform (lambda (y) (/ x y))
                            average-damp
                            1.0))

;(define (sqrt x)
;  (fixed-point-of-transform (lambda (y) (- (square y) x))
;                            newton-transform
;                            1.0))


;;;EXERCISE 1.40

(define (cubic a b c)
  (lambda (x)
    (+ (cube x)
       (* a (square x))
       (* b x)
       c)))
;Test
;(newtons-method (cubic 3 2 1) 1)

;;;EXERCISE 1.41

;(define (double f)
;    (lambda (x)
;        (f (f x))))

(define ((double f) x)
  (f (f x)))

(((double (double double)) inc) 5)

;;;EXERCISE 1.42
;;: ((compose square inc) 6)

(define (compose f g)
  (lambda (x)
    (f (g x))))

;;;EXERCISE 1.43

;;: ((repeated square 2) 5)

;(define (repeated f n)
;    (if (= n 1)
;        f
;        (lambda (x)
;            (let ((fs (repeated f (- n 1))))
;                (f (fs x))))))

; (define (repeated f n)
;    (if (= n 1)
;        f
;        (lambda (x)
;            (f ((repeated f (- n 1)) x)))))

;(define (repeated f n)
;    (define (iter i repeated-f)
;        (if (= i 1)
;            repeated-f
;            (iter (- i 1)
;                  (lambda (x)
;                      (f (repeated-f x))))))
;    (iter n f))

;(define (repeated f n)
;    (if (= n 1)
;        f
;        (compose f
;                 (repeated f (- n 1)))))

(define (repeated f n)
  (define (iter i repeated-f)
    (if (= i 1)
        repeated-f
        (iter (- i 1)
              (compose f repeated-f))))
  (iter n f))


;;;EXERCISE 1.44

(define (smooth f)
  (lambda (x)
    (/ (+ (f (- x dx))
          (f x)
          (f (+ x dx)))
       3)))

;(define (smooth-n-times f n)
;    (if (= n 0)
;        f
;        (smooth (smooth-n-times f (- n 1)))))

;(define (smooth-n-times f n)
;    (define (iter i smoothed-f)
;        (if (= i 0)
;            smoothed-f
;            (iter (- i 1)
;                  (smooth smoothed-f))))
;    (iter n f))

(define (smooth-n-times f n)
  (let ((n-times-smooth (repeated smooth n)))
    (n-times-smooth f)))


;;;EXERCISE 1.45

(define (expt base n)
  (if (= n 0)
      1
      ((repeated (lambda (x) (* base x)) n) 1)))

(define (average-damp-n-times f n)
  ((repeated average-damp n) f))

(define (damped-nth-root n damp-times)
  (lambda (x)
    (fixed-point 
     (average-damp-n-times 
      (lambda (y) 
        (/ x (expt y (- n 1)))) 
      damp-times)
     1.0)))

(define (lg n)
  (if (= n 1)
      1
      (cond ((> (/ n 2) 1)
             (+ 1 (lg (/ n 2))))
            ((< (/ n 2) 1)
             0)
            (else
             1))
      ))

(define (nth-root n)
  (damped-nth-root n (lg n)))

;;;EXERCISE 1.45

(define (iterative-improve close-enough? improve)
  (lambda (first-guess)
    (define (try guess)
      (let ((next (improve guess)))
        (if (close-enough? guess next)
            next
            (try next))))
    (try first-guess)))

;(define (fixed-point f first-guess)
;    (define tolerance 0.00001)
;    (define (close-enough? v1 v2)
;        (< (abs (- v1 v2)) tolerance))
;    (define (improve guess)
;        (f guess))
;    ((iterative-improve close-enough? improve) first-guess))
;
;(define (sqrt x)
;    (define dx 0.00001)
;    (define (close-enough? v1 v2)
;        (< (abs (- v1 v2)) dx))
;    (define (improve guess)
;        (average guess (/ x guess)))
;    (define (average x y)
;        (/ (+ x y) 2))
;    ((iterative-improve close-enough? improve) 1.0))