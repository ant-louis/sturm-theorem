#lang racket

(provide p%q)
(provide sturm-chain)
(provide count-roots)
(provide find-roots)

; All polynomes are represented by the list of size (d+1) of their coefficients
; in increasing order of exponent, where d is the highest exponent whose
; coefficient is non-zero.
; For instance, the polynme P(x) = 2 + 3x + x^3  is represented by (2 3 0 1).


;******************************MODULO_OF_TWO_POLYNOMES*******************************

; If p and q are respectively the representation of P(x) and Q(x),
; (p%q p q) returns the representation of R(x), the remainer of the division
; of P(x) by Q(x).
(define p%q
  (lambda (p q)
    (let ([revr (p%q-aux (reverse p) (reverse q))])
              (reverse revr))))

; If p and q are respectively the reversed representation of P(x) and Q(x) (means the 
; representation by decreasing order of exponent), (p%q-aux p q) returns the inversed
; representation of R(x), the remainer of the division of P(x) by Q(x).
(define p%q-aux
  (lambda (p q)
    (cond ((= (length p) 1) p)
          ((= (car p) 0) (p%q-aux (cdr p) q))
          ((< (length p) (length q)) p)
          ((let* ([div (/ (car p) (car q))]
                  [tmp (map(lambda(x) (* div x)) q)]
                  [new (sub-pol p tmp)])
            (p%q-aux (cdr new) q))))))

; If p and q are respectively the representation of P(x) and Q(x),
; (sub-pol p q) returns the representation of R(x), the result of the substraction
; of P(x) by Q(x).
(define sub-pol
  (lambda (p q)
    (cond ((null? q) p)
          ((null? p) q)
          (else (cons (- (car p) (car q)) (sub-pol (cdr p) (cdr q)))))))



;************************************STURM_CHAIN***************************************

; If p is the representation of P(x), (sturm-chain p) returns the Sturm chain of p by
; decreasing order of powers.
(define sturm-chain
  (lambda (p)
    (cond ((<= (length p) 1) (list p))
          ((= (length p) 2) (list p (deriv p)))
          (else (sturm-chain-aux (list p (deriv p)))))))

; If ls is a list of n polynomes (n >= 2), (sturm-chain-aux ls) returns the
; concatenation of ls and the sequence of the k non-null remainers of the division
; where ls_k = ls_k-2 % ls_k-1.
(define sturm-chain-aux
  (lambda (ls)
    (let* ([n (length ls)]
           [j-1 (list-ref ls (- n 1))]
           [j-2 (list-ref ls (- n 2))]
           [j  (map - (p%q j-2 j-1))])
      (if (<= (length j) 1) (append ls (list j))
          (sturm-chain-aux (append ls (list j)))))))

; If p is the representation of P(x), (deriv p) returns P'(x), the derived polynome.
(define deriv
  (lambda (p)
    (deriv-aux (cdr p) 1)))

; If p is the representation of P(x) and n the degree of the first term of that polynome (n > 0),
; (deriv-aux p n) returns P'(x), the derived polynome.
(define deriv-aux
  (lambda (p n)
         (if (null? p) '()
             (cons (* (car p) n) (deriv-aux (cdr p) (+ n 1))))))
    

  
;*************************************NB_OF_ROOTS**************************************

; If p is a polynome and both a and b are numbers such that a < b,
; (count-roots p a b) returns the number of roots of p on ]a b].
(define count-roots
  (lambda (p a b)
    (let* ([chain (sturm-chain p)]
           [eval_a (eval-seq chain a)]
           [eval_b (eval-seq chain b)]
           [nb_a (sign-change eval_a)]
           [nb_b (sign-change eval_b)])
      (- nb_a nb_b))))

; If ls is a list of numbers, (sign-change ls) returns the number of sign changes
; occuring in the list beginning from the first element of ls.
(define sign-change
  (lambda (ls)
    (sign-change-aux (cdr ls) (car ls))))

; If ls is a list of numbers and a is a number, (sign-change-aux ls a) returns
; the number of sign changes in the sequence from a as beginning point.
(define sign-change-aux
  (lambda (ls a)
    (cond ((null? ls) 0)
          ((< (* a (car ls)) 0) (+ 1 (sign-change-aux (cdr ls) (car ls))))
          (else (sign-change-aux (cdr ls) (car ls))))))
           
; If ls is the representation of the Sturm chain of a polynome P(x) and a is a number, 
; (eval-seq ls a) returns the sequence of the evaluations in x=a for each term of the
; given Sturm chain.      
(define eval-seq
  (lambda (ls a)
    (map (lambda (x) (eval-pol x a 0)) ls)))

; If p is the representation of a polynome P(x), a is a number and n the minimal degree of P(x)
; (n >= 0), (eval-pol p a n) returns the evaluation of p in x=a.
(define eval-pol
  (lambda (p a n)
    (cond ((null? p) 0)
          ((or (= a +inf.0) (= a -inf.0)) (* (car (reverse p)) (expt a (- (length p) 1))))
          (else (+ (* (car p) (expt a n)) (eval-pol (cdr p) a (+ n 1)))))))



;***********************************FIND_ROOTS(BONUS)**********************************
    
; If p is a polynome, both a and b are numbers (such that a < b) and eps
; is a positive real, (find-roots p a b eps) returns the ordered list
; of roots of p on the ]a, b] interval with precision eps.
(define find-roots
  (lambda (p a b eps)
    (find-roots-aux p a b eps '())))
                       
; If p is a polynome, both a and b are numbers (such that a < b), eps
; is a positive real and ls a list, (bisection p a b eps) returns the 
; concatenation of ls and the ordered list of roots of p on the ]a, b]
; interval with precision eps.
(define find-roots-aux
  (lambda (p a b eps ls)
    (let ([m (/ (+ a b) 2)]
          [nb (count-roots p a b)])
         (cond ((zero? nb) ls)
               ((and (= nb 1) (< (/ (- b a) 2) eps)) (cons m ls))
               (else (find-roots-aux p a m eps (find-roots-aux p m b eps ls)))))))
                    