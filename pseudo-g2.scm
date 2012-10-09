;;; Pseudo G2 Test
;;; by Peter Salvi, 2012

;;; Sensible configurations
;;; [L: linear, P: parabolic/cubic-linear, C: parabolic/common, G: gamma, Q: quintic]
;;;
;;; [1] L..
;;; Two linear ribbons blended by a cubic Hermite function.
;;; [2] LG.
;;; Same as above, but reparameterized such that it interpolates
;;; the cubic Hermite curve.
;;; [3] P.Q
;;; Two parabolic ribbons blended by a quintic Hermite function.
;;; The ribbons' radii of curvature are the same as in [1].
;;; [4] PGQ
;;; Same as above, but the radii of curvature are the same as in [2].
;;; [5] C..
;;; Two parabolic ribbons blended by a cubic Hermite function.
;;; The ribbons' radii of curvature are set in such a way,
;;; that they interpolate the value set on the slider.
;;; [6] C.Q
;;; Two parabolic ribbons blended by a quintic Hermite function.
;;; The ribbons' radii of curvature are set on the slider.
;;; [7] CGQ
;;; Same as above, but reparameterized such that it interpolates
;;; the quintic Hermite curve.

#lang racket

;;; Libraries
(require racket/gui)

;;; Parameters
(define alpha 1/4)         ; placement of the curved ribbon's first control point
(define point-radius 4)
(define line-width 2)
(define resolution 100)
(define curvature-resolution 60)
(define curvature-scaling 1/3)
(define curvature-fit-points 5) ; on either side

;;; Default placements
(define p0 '(230.0 480.0))
(define p1 '(480.0 190.0))
(define p2 '(730.0 480.0))
(define d0 '(0.0 -300.0))
(define d1a '(-300.0 0.0))
(define d1b '(300.0 0.0))
(define d2 '(0.0 -300.0))

;;; Options
(define points? #t)
(define circle? #f)
(define linear? #t)
(define parabolic? #f)
(define curved? #f)
(define control-points? #f)
(define curvature? #f)
(define curvature 300)
(define gamma? #t)
(define quintic? #f)
(define simple-mean? #t)

;;; Variables
(define dragged #f)
(define deviation #f)

;;; Basic Maths
(define (binomial n k)
  (if (= k 0)
      1
      (* (/ n k) (binomial (- n 1) (- k 1)))))
(define (v+ . args) (apply map + args))
(define (v- . args) (apply map - args))
(define (v* u . args) (map (lambda (x) (apply * x args)) u))
(define (vlength u) (sqrt (apply + (map (lambda (x) (* x x)) u))))
(define (vnormalize u) (v* u (/ (vlength u))))
(define (vperp u) (list (second u) (- (first u))))
(define (point-distance p q) (vlength (v- q p)))
(define (scalar-product u v) (apply + (map * u v)))
(define (cross-product-length u v) (abs (- (* (first u) (second v)) (* (second u) (first v)))))
(define (rotate-vector d p)
  (let* ([u (vnormalize d)]
         [v (vperp u)])
    (list (scalar-product p u) (scalar-product p v))))
(define (rotate-to o d p)
  (rotate-vector d (v- p o)))

;;; Curves

;;; gamma = H1/H0, where
;;;   H0 = 2u3-3u2+1 or -6u5+15u4-10u3+1 [see (hermite 0)]
;;; and
;;;   H1 = u3-2u2+u or -3u5+8u4-6u3+u
(define (gamma u)
  (if gamma?
      (if quintic?
          (/ (* u (+ (* 3 u) 1))
             (+ (* 6 u u) (* 3 u) 1))
          (/ u (+ (* 2 u) 1)))
      u))

(define (gamma-derivative k u)
  (if gamma?
      (if quintic?
          (let ([denom (+ (* 6 u u) (* 3 u) 1)])
            (if (= k 1)
                (/ (+ (* 3 u u) (* 6 u) 1) denom denom)
                (/ (* -36 u (+ (* u u) (* 3 u) 1)) denom denom denom)))
          (let ([denom (+ (* 2 u) 1)])
            (if (= k 1)
                (/ 1 denom denom)
                (/ -4 denom denom denom))))
      (if (= k 1) 1 0)))

;;; delta = H2/H0, where
;;;   H2 = -1/2u5+3/2u4-3/2u3+1/2u2
(define (delta u)
  (if (and gamma? quintic?)
      (/ (* u u)
         (* 2 (+ (* 6 u u) (* 3 u) 1)))
      (* u u 1/2)))

(define (delta-derivative k u)
  (if (and gamma? quintic?)
      (let ([denom (+ (* 6 u u) (* 3 u) 1)])
            (if (= k 1)
                (/ (* u (+ (* 3 u) 2)) 2 denom denom)
                (/ (+ 1 (* -18 u u u) (* -18 u u)) denom denom denom)))
      (if (= k 1) u 1)))

(define (linear-ribbon p d)
  (lambda (u)
    (v+ p (v* d (gamma u)))))

(define (linear-ribbon-derivative d)
  (lambda (k u) (v* d (gamma-derivative k u))))

(define (curved-ribbon-cp-% p d c)
  (let ([n (vnormalize (vperp d))]
        [q (v* (v+ p0 p2) 1/2)])
    (when (< (scalar-product (v- q p1) n) 0)
      (set! n (v* n -1)))
    (let ([correction
           (if (or quintic? (not (or (equal? d d1a) (equal? d d1b))))
               0
               (let* ([xp (rotate-to p1 n (if (equal? d d1a) p0 p2))]
                      [xt (rotate-vector n (v* (if (equal? d d1a) d0 d2) 2 alpha))]
                      [xq (v* xt (- (/ alpha) 2))])
                 (/ (+ (* 6 (first xp)) (* 6 (first xt)) (* 3 (first xq))) 2)))])
      (list p (v+ p (v* d alpha))
            (if (= c 0)
                (v+ p d)
                (v+ p d (v* n (- (* (scalar-product d d) 2 alpha alpha (/ c)) correction))))))))

(define (curved-ribbon p d c)
  (let* ([points (curved-ribbon-cp-% p d c)]
         [q0 (first points)]
         [q1 (v* (v- (second points) (first points)) 2)]
         [q2 (v* (v+ (first points) (v* (second points) -2) (third points)) 2)])
    (lambda (u)
      (v+ q0 (v* q1 (gamma u)) (v* q2 (delta u))))))

(define (curved-ribbon-derivative p d c)
  (let* ([points (curved-ribbon-cp-% p d c)]
         [q1 (v* (v- (second points) (first points)) 2)]
         [q2 (v* (v+ (first points) (v* (second points) -2) (third points)) 2)])
    (lambda (k u)
      (v+ (v* q1 (gamma-derivative k u)) (v* q2 (delta-derivative k u))))))

(define (curved-ribbon-points p d c)
  (let ([r (curved-ribbon p d c)])
    (for/list ([i (in-range resolution)])
      (let* ([u (/ i (- resolution 1) 3)]
             [p (r u)])
        (make-object point% (first p) (second p))))))

(define (curved-ribbon-cp p d c)
  (map (lambda (p) (make-object point% (first p) (second p)))
       (curved-ribbon-cp-% p d c)))

(define (hermite k)
  (if quintic?
      (case k
        [(0) (lambda (u) (+ (* -6 (expt u 5)) (* 15 (expt u 4)) (* -10 u u u) 1))]
        [(1) (lambda (u) (+ (* -30 (expt u 4)) (* 60 u u u) (* -30 u u)))]
        [(2) (lambda (u) (+ (* -120 u u u) (* 180 u u) (* -60 u)))])
      (case k
        [(0) (lambda (u) (+ (* 2 u u u) (* -3 u u) 1))]
        [(1) (lambda (u) (+ (* 6 u u) (* -6 u)))]
        [(2) (lambda (u) (+ (* 12 u) -6))])))

(define (curve-points ribbon1 ribbon2)
  (let ([blend (hermite 0)])
    (for/list ([i (in-range resolution)])
      (let* ([u (/ i (- resolution 1))]
             [p (v+ (v* (ribbon1 u) (blend u))
                    (v* (ribbon2 (- 1 u)) (blend (- 1 u))))])
        (make-object point% (first p) (second p))))))

(define (curve-derivative r1 r1d r2 r2d k)
  (let ([blend (hermite 0)]
        [blend-d1 (hermite 1)]
        [blend-d2 (hermite 2)])
    (if (= k 1)
        (lambda (u)
          (v+ (v* (r1d 1 u) (blend u))
              (v* (r1 u) (blend-d1 u))
              (v* (r2d 1 (- 1 u)) (- (blend (- 1 u))))
              (v* (r2 (- 1 u)) (- (blend-d1 (- 1 u))))))
        (lambda (u)
          (v+ (v* (r1d 2 u) (blend u))
              (v* (r1d 1 u) 2 (blend-d1 u))
              (v* (r1 u) (blend-d2 u))
              (v* (r2d 2 (- 1 u)) (blend (- 1 u)))
              (v* (r2d 1 (- 1 u)) 2 (blend-d1 (- 1 u)))
              (v* (r2 (- 1 u)) (blend-d2 (- 1 u))))))))

(define (curve-curvature r1 r1d r2 r2d)
  (lambda (u)
    (let* ([d1 ((curve-derivative r1 r1d r2 r2d 1) u)]
           [d2 ((curve-derivative r1 r1d r2 r2d 2) u)]
           [l (vlength d1)]
           [d1d2 (cross-product-length d1 d2)])
      (if (not (= d1d2 0))
          (/ (* l l l) d1d2)
          0))))

(define (show-segment dc a b)
  (send dc draw-lines
        (list (make-object point% (first a) (second a))
              (make-object point% (first b) (second b)))))

(define (show-curvature-lines dc r1 r1d r2 r2d)
  (let ([der (curve-derivative r1 r1d r2 r2d 1)]
        [curv (curve-curvature r1 r1d r2 r2d)]
        [blend (hermite 0)])
    (for/list ([i (in-range curvature-resolution)])
      (let* ([u (/ i (- curvature-resolution 1))]
             [p (v+ (v* (r1 u) (blend u))
                    (v* (r2 (- 1 u)) (blend (- 1 u))))]
             [n (vnormalize (vperp (der u)))]
             [c (curv u)])
        (show-segment dc p (v+ p (v* n curvature-scaling c)))))))

(define (curvature-center)
  (let ([n (vnormalize (vperp d1a))]
        [q (v* (v+ p0 p2) 1/2)])
    (if (< (scalar-product (v- q p1) n) 0)
        (v+ p1 (v* n -1 curvature))
        (v+ p1 (v* n curvature)))))

(define (gather-points)
  (let ([gather (lambda (q1 v1 q2 v2)
                  (let ([blend (hermite 0)]
                        [r1 (linear-ribbon q1 v1)]
                        [r2 (linear-ribbon q2 v2)])
                    (for/list ([i (in-range curvature-fit-points)])
                      (let ([u (/ (+ i 1) (- curvature-resolution 1))])
                        (v+ (v* (r1 u) (blend u))
                            (v* (r2 (- 1 u)) (blend (- 1 u))))))))])
    (append (gather p1 d1a p0 d0) (gather p1 d1b p2 d2))))

(define (circle-fit)
  (let* ([points (map (lambda (p) (rotate-to p1 d1b p)) (gather-points))]
         [radii (map (lambda (p)
                       (let ([x (first p)] [y (second p)])
                         (/ (+ (* x x) (* y y)) (* 2 y))))
                     points)])
    (abs (/ (foldl + 0 radii) (length radii)))))

(define (curvatures-left-right)
  (let ([q quintic?])
    (set! quintic? #f)
    (let* ([left ((curve-curvature (linear-ribbon p0 d0) (linear-ribbon-derivative d0)
                                   (linear-ribbon p1 d1a) (linear-ribbon-derivative d1a))
                  1)]
           [right ((curve-curvature (linear-ribbon p1 d1b) (linear-ribbon-derivative d1b)
                                    (linear-ribbon p2 d2) (linear-ribbon-derivative d2))
                   0)]
           [result (list left right)])
      (set! quintic? q)
      result)))

(define (circle-mean)
  (let ([lr (curvatures-left-right)])
    (/ (+ (first lr) (second lr)) 2)))

(define (best-curvature)
  (if (or quintic? simple-mean?)
      (circle-mean)
      (circle-fit)))

(define (curvature-deviation)
  (let ([left ((curve-curvature (curved-ribbon p0 d0 0) (curved-ribbon-derivative p0 d0 0)
                                (curved-ribbon p1 d1a curvature) (curved-ribbon-derivative p1 d1a curvature)) 1)]
        [right ((curve-curvature (curved-ribbon p1 d1b curvature) (curved-ribbon-derivative p1 d1b curvature)
                                 (curved-ribbon p2 d2 0) (curved-ribbon-derivative p2 d2 0)) 0)])
    (if (= curvature 0)
        (list 0 0)
        (map (lambda (x) (inexact->exact (round (* (/ x curvature) 100))))
             (list (- curvature left) (- curvature right))))))

;;; Graphics

(define (draw-point dc p)
  (send dc draw-ellipse
        (- (first p) point-radius) (- (second p) point-radius)
        (* point-radius 2) (* point-radius 2)))

(define (draw canvas dc)
  (when linear?
    (send dc set-pen "GREEN" line-width 'solid)
    (send dc draw-lines (curve-points (linear-ribbon p0 d0) (linear-ribbon p1 d1a)))
    (send dc draw-lines (curve-points (linear-ribbon p1 d1b) (linear-ribbon p2 d2)))
    (when (and curvature? (not dragged))
      (show-curvature-lines dc (linear-ribbon p0 d0) (linear-ribbon-derivative d0)
                            (linear-ribbon p1 d1a) (linear-ribbon-derivative d1a))
      (show-curvature-lines dc (linear-ribbon p1 d1b) (linear-ribbon-derivative d1b)
                            (linear-ribbon p2 d2) (linear-ribbon-derivative d2))))
  (when parabolic?
    (send dc set-pen "PINK" line-width 'solid)
    (let* ([lr (curvatures-left-right)]
           [left (first lr)] [right (second lr)])
      (send dc draw-lines (curve-points (curved-ribbon p0 d0 0) (curved-ribbon p1 d1a left)))
      (send dc draw-lines (curve-points (curved-ribbon p1 d1b right) (curved-ribbon p2 d2 0)))
      (when (and curvature? (not dragged))
        (show-curvature-lines dc (curved-ribbon p0 d0 0) (curved-ribbon-derivative p0 d0 0)
                              (curved-ribbon p1 d1a left) (curved-ribbon-derivative p1 d1a left))
        (show-curvature-lines dc (curved-ribbon p1 d1b right) (curved-ribbon-derivative p1 d1b right)
                              (curved-ribbon p2 d2 0) (curved-ribbon-derivative p2 d2 0)))))
  (when curved?
    (send dc set-pen "LIGHT BLUE" line-width 'solid)
    (send dc draw-lines (curve-points (curved-ribbon p0 d0 0) (curved-ribbon p1 d1a curvature)))
    (send dc draw-lines (curve-points (curved-ribbon p1 d1b curvature) (curved-ribbon p2 d2 0)))
    (when (and curvature? (not dragged))
      (show-curvature-lines dc (curved-ribbon p0 d0 0) (curved-ribbon-derivative p0 d0 0)
                            (curved-ribbon p1 d1a curvature) (curved-ribbon-derivative p1 d1a curvature))
      (show-curvature-lines dc (curved-ribbon p1 d1b curvature) (curved-ribbon-derivative p1 d1b curvature)
                            (curved-ribbon p2 d2 0) (curved-ribbon-derivative p2 d2 0))))
  (when (and circle? (not (= curvature 0)))
    (send dc set-pen "BLACK" line-width 'short-dash)
    (send dc set-brush "BLACK" 'transparent)
    (let* ([2r (* curvature 2)]
           [center (curvature-center)]
           [left (- (first center) curvature)]
           [top (- (second center) curvature)])
      (send dc draw-ellipse left top 2r 2r)))
  (when control-points?
    (send dc set-pen "RED" line-width 'solid)
    (let ([g gamma?])
      (set! gamma? #f)
      (send dc draw-lines (curved-ribbon-cp p1 d1a curvature))
      (send dc draw-lines (curved-ribbon-points p1 d1a curvature))
      (send dc draw-lines (curved-ribbon-cp p1 d1b curvature))
      (send dc draw-lines (curved-ribbon-points p1 d1b curvature))
      (set! gamma? g)))
  (when points?
    (let ([p0+d0 (v+ p0 (v* d0 1/3))]
          [p1+d1a (v+ p1 (v* d1a 1/3))]
          [p1+d1b (v+ p1 (v* d1b 1/3))]
          [p2+d2 (v+ p2 (v* d2 1/3))])
      (send dc set-brush "BLACK" 'solid)
      (send dc set-pen "BLACK" line-width 'solid)
      (for-each (lambda (p) (draw-point dc p)) (list p0 p1 p2))
      (send dc set-brush "BLUE" 'solid)
      (send dc set-pen "BLUE" line-width 'solid)
      (for-each (lambda (p) (draw-point dc p)) (list p0+d0 p1+d1a p1+d1b p2+d2))
      (send dc draw-line (first p0) (second p0) (first p0+d0) (second p0+d0))
      (send dc draw-line (first p1) (second p1) (first p1+d1a) (second p1+d1a))
      (send dc draw-line (first p1) (second p1) (first p1+d1b) (second p1+d1b))
      (send dc draw-line (first p2) (second p2) (first p2+d2) (second p2+d2))))
  (let ([dev (curvature-deviation)])
    (send deviation set-label (format "Curvature deviation: [left] ~a% [right] ~a%"
                                      (first dev) (second dev)))))

;;; GUI

(define (handle-mouse-movement event)
  (if dragged
      (let ([p (list (send event get-x) (send event get-y))])
        (case dragged
          [(0) (set! p0 p)]
          [(1) (set! d0 (v* (v- p p0) 3))]
          [(2) (set! p1 p)]
          [(3) (set! d1a (v* (v- p p1) 3))
               (set! d1b (v* (vnormalize d1a) (vlength d1b) -1))]
          [(4) (set! d1b (v* (v- p p1) 3))
               (set! d1a (v* (vnormalize d1b) (vlength d1a) -1))]
          [(5) (set! p2 p)]
          [(6) (set! d2 (v* (v- p p2) 3))])
        #t)
    #f))

(define (handle-mouse-down event)
  (if dragged
      (handle-mouse-up event)
      (let* ([p (list (send event get-x) (send event get-y))]
             [p0+d0 (v+ p0 (v* d0 1/3))]
             [p1+d1a (v+ p1 (v* d1a 1/3))]
             [p1+d1b (v+ p1 (v* d1b 1/3))]
             [p2+d2 (v+ p2 (v* d2 1/3))]
             [points (list p0 p0+d0 p1 p1+d1a p1+d1b p2 p2+d2)])
        (for ([i (in-range 7)]
              [pi points])
          (when (< (point-distance p pi) point-radius)
            (set! dragged i))))))

(define (handle-mouse-up event)
  (set! dragged #f)
  #t)

(define my-canvas%
  (class canvas%
    (inherit refresh)
    (define/override (on-event event)
      (when (case (send event get-event-type)
              [(motion) (handle-mouse-movement event)]
              [(left-down) (handle-mouse-down event)]
              [(left-up) (handle-mouse-up event)])
        (refresh)))
    (super-new)))

(let* ([frame (new frame% [label "Linear vs. Curved Ribbons"])]
       [vbox (new vertical-pane% [parent frame])]
       [canvas (new my-canvas% [parent vbox]
                    [min-width 800]
                    [min-height 600]
                    [paint-callback draw])]
       [hbox1 (new horizontal-pane% [parent vbox])]
       [hbox2 (new horizontal-pane% [parent vbox])]
       [hbox3 (new horizontal-pane% [parent vbox])])
  (let-syntax ([add-option (syntax-rules ()
			     [(_ var a-label box)
			      (new check-box% [parent box]
                                   [label a-label] [value var]
				   [callback (lambda (c e)
					       (set! var (not var))
					       (send canvas refresh))])])])
    (add-option points? "Points" hbox1)
    (add-option linear? "Curve (linear)" hbox1)
    (add-option parabolic? "Curve (parabolic - curv. from cubic linear)" hbox1)
    (add-option curved? "Curve (parabolic - common curvature)" hbox1)
    (add-option control-points? "Ribbon (parabolic - common curvature)" hbox1)
    (add-option circle? "Show circle" hbox2)
    (add-option curvature? "Show curvature" hbox2)
    (add-option gamma? "Use gamma function" hbox3)
    (add-option quintic? "Use quintic weights" hbox3)
    (add-option simple-mean? "Use simple mean [always on in quintic mode]" hbox3)
    (set! deviation (new message% [label ""] [parent hbox3]
                         [stretchable-width #t])))
  (let ([slider
         (new slider% [label "Curvature radius:"] [parent hbox2]
              [style '(horizontal plain)] [init-value curvature]
              [min-value 0] [max-value 1000]
              [callback (lambda (b e)
                          (set! curvature (send b get-value))
                          (send canvas refresh))])])
    (new button% [parent hbox2] [label "Best fit curvature"] [stretchable-width #t]
         [callback (lambda (b e)
                     (set! curvature (best-curvature))
                     (send canvas refresh)
                     (send slider set-value (inexact->exact (round curvature))))]))
  (send frame show #t))
