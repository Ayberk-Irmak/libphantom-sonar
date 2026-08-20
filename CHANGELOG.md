# Changelog

Every entry records what was measured, and — where it applies — what was
found to be wrong. Several releases exist because a test disagreed with the
code and the code turned out to be right; several others exist because a
measurement was being compared against itself. Both are kept.

Written newest first. The full reasoning for any number here is in
[`docs/math_spec.md`](docs/math_spec.md) and
[`docs/validation.md`](docs/validation.md).

## 0.16.0 — turning a recording into a measurement

v0.15 built an air bench and could not run it: this machine has no
working acoustic path. How that was FOUND turned out to matter more than
the fact of it, and this release is the machinery that came out of it.

It was not found by looking at the recording. The recording looked fine
-- 13000 distinct sample values, a healthy RMS, a matched filter peak
44 dB above the background. Every one of those is what a working
microphone gives, and every one was produced by electrical noise and a
startup transient. What settled it was a CONTROLLED comparison: record
while playing, record in silence, difference the in-band energy. The
answer came out NEGATIVE -- playing made it quieter, which no working
channel can do.

QUALIFICATION, BEFORE ANY NUMBER IS TAKEN.

    channel               excess   clipped   verdict
    carrying the probe   +18.92dB      0%    usable
    carrying nothing      -0.27dB      0%    unusable
    carrying it, loudly  +31.26dB   18.24%   unusable

The third row is the one that earns the test. It passes the "is the
signal there" check by a wide margin and is rejected anyway, because a
clipping recording reports levels that are fiction. A qualifier that only
asked whether the probe arrived would wave it through.

Band energy is measured through a Hann window, because without one the
spectral leakage from a loud out-of-band transient lands squarely in the
band being measured -- which is precisely how a startup click passes for
signal.

TWO DISTANCES CANCEL THE SOUND CARD.

A delay measured over a desk is almost entirely buffering: 10-50 ms of it
against 1.5 ms of flight time over half a metre.

    naive single-shot d1/t1   12.1 m/s   -- out by 28x
    two-distance solve        343.37 m/s, latency 32.00 ms

Both recovered exactly. The solver refuses distances too close together
to resolve, rather than dividing by a difference that is mostly timing
noise.

IMPULSE RESPONSE BY REGULARISED DECONVOLUTION.

Three known arrivals at 200/320/512 samples with gains 1.0/0.5/0.25,
recovered at 200.0/320.0/512.0 and -6.03/-12.06 dB against truths of
-6.02/-12.04.

The regularisation is not a refinement. A probe has almost no energy
outside its own band, so dividing by |S|^2 there amplifies whatever noise
is present without bound -- and the result still looks like an impulse
response. The floor is a fraction of the probe's own mean power rather
than an absolute constant, so the estimate does not depend on the
recording's gain.

Also fixed from v0.15: the arrival picker's separation guard was 20 ms,
which is seven metres of extra path and silently discarded every echo a
room produces -- a direct path and its first reflection are typically
2-4 ms apart. It is now a few times the probe's own 1/B resolution.

STILL NOT MEASURED, and every document says so. No number in this project
has been taken from a real transducer. This is the apparatus, verified
against channels whose answers are known by construction; the
qualification step is what will decide whether a real setup can produce a
number worth having.

189 cases / 95846 checks, 269 C checks, 8 Rust tests. Clean under gcc and
clang in both precisions, under ASan+UBSan, zero allocation, and ARM
still reports data == 0 and bss == 0. Bellhop cross-validation passes.

## 0.15.0 — acoustics in air, and a bench you can actually run

Hardware in the loop, started in the medium anyone can test in. A
hydrophone, a projector and an amplifier cost several hundred euros and
need water to put them in; a speaker and a microphone are already on the
desk. The signal processing does not care which medium it is in -- a
matched filter, a Doppler bank, a CFAR detector and a tracker are the
same code at 340 m/s as at 1500 -- so an air bench exercises the whole
chain before any wet hardware is bought.

ATMOSPHERIC ABSORPTION, AGAINST THE STANDARD ITSELF.

ISO 9613-1:1993 equations (3) to (5), verified against Table 1 of the
same document: 105 published values at 10 C, worst relative error 0.380%,
and that worst case is at 315 Hz where the table prints 1.30 against a
computed 1.305 -- the table's own three-figure rounding rather than a
disagreement.

The transcription was checked the way v0.11's UNESCO table was: the
equations entered separately from clause 6.2, so a slipped digit appears
as one large outlier rather than a uniform near-miss. The PDF's text
layer is OCR-damaged in places, so only unambiguous cells were taken.

Note 5 of the standard, quantified. It records that Table 1 was computed
at the EXACT one-third-octave midband frequencies 1000*10^(k/10), not the
preferred nominal values printed in its own row headings. Using the
nominal ones makes the worst error 1.733% instead of 0.380% -- 4.6x
worse, small, systematic, and exactly the kind of thing that gets blamed
on an implementation.

HUMIDITY IS NOT A CORRECTION, IT IS THE EFFECT.

At 10 C and 4 kHz, going from 10% to 20% relative humidity MULTIPLIES the
loss by 1.60. And the sign reverses with frequency: at 500 Hz damp air
absorbs less, at 8 kHz more. A model that only knew "damp air absorbs
more" would be wrong half the time. Nothing in seawater behaves like
this, and it is why a bench measurement in air that does not record the
humidity means nothing.

AIR IS THE HARSHER DOPPLER ENVIRONMENT, BY 4.37x.

v/c is 2.91e-3 in air against 6.67e-4 in water. A 511-chip spreading code
slips 1.49 chips per bit at walking pace, where in water it slips 0.34.
That makes an air bench a STRONGER test of v0.13's comm module than a
water one, not a weaker substitute.

Also: sound speed (343.37 m/s at 20 C dry, and the header says it is not
Cramer 1993 and why that is fine here), and impedance -- 412 rayl against
seawater's 1.5e6, a factor of 3637, which is why an air bench says
nothing quantitative about target strength however well it exercises the
processing.

AND THE HONEST PART.

Real hardware-in-the-loop was attempted on this machine and does not
work. The audio devices enumerate -- an ALC256 with analogue playback and
capture -- but there is no working acoustic path. A controlled test
settles it: recording while playing five 50 ms chirps gave 7.85 dB LESS
energy in the 2-8 kHz band than recording in silence, and capture is
dominated by a large startup transient. No speaker-to-microphone coupling
could be demonstrated, so no hardware measurement is claimed anywhere.

tools/air_bench.py is written and self-verified instead: its selftest
runs the full analysis against a simulated channel and recovers a direct
path at 1.500 ms (truth 1.500) and an echo at 4.00 ms, 8 dB down. Two
traps are documented in the tool rather than left to be found -- it
measures timing and detection but not level, since consumer transducers
are uncalibrated; and the measured delay includes the sound card's own
10-50 ms latency, which dwarfs the 1.5 ms sound takes to cross half a
metre, so only a two-distance measurement removes it.

Two smaller fixes found on the way: the bench's peak picker used a 20 ms
separation and silently discarded every echo arriving within 7 metres of
the direct path, which on a desk is all of them; and its sound-speed
helper defaulted to 50% humidity where the C++ defaults to dry, so two
things documented as "kept in step" disagreed by 0.6 m/s. Also
build.rs now watches the archive, since cargo was happily relinking
against a stale library and testing the old binary.

186 cases / 95825 checks, 269 C checks, 8 Rust tests. Clean under gcc and
clang in both precisions, under ASan+UBSan, zero allocation. Bellhop
cross-validation still passes.

## 0.14.0 — bindings and portability

No new physics. This release is about whether the physics can leave the
building, and it found one real defect and one honest gap.

A C ABI WHERE THE CALLER OWNS THE MEMORY.

Hand-written, not generated: a generator mirrors whatever the headers say
today including their mistakes, and the C++ interface is spans, templates
and RAII, none of which have a stable C representation. The library
allocates nothing and the ABI does not quietly change that by handing out
pointers a caller must free -- every stateful object is placed into
storage the caller supplies, with the library deciding only how big it is.

Tested by a C COMPILER in C11 mode with -Wstrict-prototypes -Werror, not
by C++ in C mode. A header that only ever meets a C++ compiler
accumulates C++-isms and nobody finds out until a real C caller arrives.
269 checks, round-tripping values with published answers -- the UNESCO
check value, the CRC-32 check value -- rather than merely returning them.

160 kB OF .BSS THAT CROSS-COMPILING FOUND.

The first C ABI copied traced rays through a static 8192-point scratch
buffer, on the reasoning that it "must not assume ph_ray_point and
RayPoint share a layout". Building for ARM and reporting section sizes
showed what that cost: 163840 bytes of .bss, half the RAM of the
Cortex-M7 this library is meant to fit on, for a copy that does nothing.
It also made ph_trace_ray non-reentrant.

The assumption did not need avoiding, it needed checking. static_asserts
on size, alignment and every member offset now break the BUILD if it ever
stops holding, instead of the ABI silently scrambling ray paths. With the
tracker's measurement conversion moved from a static to a stack buffer:

    text 61216 -> 61160
    data  1536 -> 0
    bss 163840 -> 0

The library now carries no mutable static storage at all, which is both
the RAM and the reentrancy. CI asserts data == 0 && bss == 0.

CROSS-COMPILATION, VERIFIED AND UNVERIFIED.

Verified: arm-linux-gnueabihf 32-bit with Real=float, and
riscv64-linux-gnu 64-bit with Real=double. Two architectures, two word
sizes, both precisions -- the combination matters because the bugs this
catches only appear when size_t, pointer width and Real change together.
Both clean under -Werror with zero allocator symbols.

NOT verified: the bare-metal Cortex-M7 build. cmake/cortex-m7.cmake is
provided and correct as far as it goes, but this environment's
arm-none-eabi toolchain ships no C++ standard library, so it could not be
run. It is labelled untested rather than claimed, and the ARM figures
above come from a hosted target with the same instruction set.

-fno-exceptions -fno-rtti builds clean and references no __cxa_throw,
_Unwind_* or __gxx_personality. A no-op that is worth having as a build:
if the library ever starts throwing, that is a link error in CI rather
than an unwinder in a firmware image.

RUST.

phantom-sonar-sys with hand-transcribed extern "C" declarations, and a
safe phantom-sonar over it. 8 tests, checking the BINDING and not the
physics -- re-asserting a ray path in Rust would prove only that FFI
copies bytes, so the tests either round-trip a published value or compare
a layout against what the library reports about itself.

The safe layer earns its place twice. A Profile owns its storage, so the
C contract's "buffer must outlive the object" becomes the borrow
checker's problem. And it can be stricter than C: ph_profile_speed_at
answers 0 m/s for a profile with fewer than two samples -- a silent
failure a caller can propagate -- where the Rust speed_at returns
Err(State). A binding test found that, not design.

Precision mismatch between library and crate is not a link error; it
reinterprets every float crossing the boundary. check_precision() exists
for that and every constructor calls it.

ALSO: tools/audit_no_alloc.py now distinguishes placement new from
allocating new. Placement new allocates nothing -- it is what makes
caller-owned storage expressible -- and banning it would ban the one
construct the C ABI is built on. The distinction is syntactic and the
tool still catches a plain `new T`.

181 C++ cases / 95806 checks, 269 C checks, 8 Rust tests. Clean under gcc
and clang in both precisions, under ASan+UBSan, zero allocation by source
audit and nm, on x86-64, ARM and RISC-V. Bellhop cross-validation still
passes.

## 0.13.0 — spread-spectrum acoustic communication

The third engine of the original specification, and the last real gap.
Same physics as the rest of the library, used the other way round: where
the ping analyser pulls a known waveform out of noise to DETECT
something, this pulls one out of noise to CARRY something.

PROCESSING GAIN, AND THE MEASUREMENT THAT WAS WRONG FIRST.

The first version of this release reported 12.17 dB of processing gain
against a theoretical 12.17 dB at 511 chips. Exact agreement, and
meaningless. It held the chip AMPLITUDE fixed while sweeping N, and
energy per bit is A^2 * N * T_chip -- so the sweep multiplied transmitted
energy by N and reported the result as gain. Both sides of the comparison
were the same tautology, which is exactly why they matched perfectly.

Corrected, and stated the narrow way that is true:

  - Against WHITE NOISE at a fixed energy per bit, spreading buys
    NOTHING. Despreading multiplies by a +/-1 sequence, which leaves
    white noise unchanged in distribution; a matched filter already
    extracts everything a known waveform in noise has to give. Measured
    flat to 0.87 dB over a 16x range of code length.
  - The real gain is BANDWIDTH EXPANSION at a fixed data rate, against
    narrowband interference. Data rate held at 100 bps and chip rate
    grown with N: +7.1 dB from 31 to 511 chips against a tone, while
    white noise stayed flat to 0.75 dB. Short of the ideal 12.2 because a
    bare correlator collects the diluted interferer from the whole band.
  - And covertness, which is the other half of why anyone spreads: the
    same power over N times the bandwidth is 10 log10(N) less spectral
    density. A transmit-side property this library does not measure, and
    says so.

A trap in this module's own API is documented at the point of use:
raising chips_per_bit at a fixed chip rate does NOT spread anything, it
trades data rate for integration time.

This is the third measurement in the project to fail by being compared
against itself, after the sound-speed equations in v0.11 and the Snell
invariant in v0.11.

THE REST.

PN sequences with maximality VERIFIED rather than assumed, degrees 5 to
15: balance exactly +1, autocorrelation exactly -1 at every shift. The
first tap table used the wrong convention for a right-shifting LFSR,
leaving bit 0 clear -- a seed of 1 fed back zero, the state collapsed on
the first step, and the generator emitted a constant with the correct
LENGTH. Only those two measurements caught it.

Doppler as a time scaling, with the ceiling it puts on gain: 10 m/s
permits 37 chips and no transmit power lifts that. HFM vs LFM preambles
measured -- at 20 m/s the HFM holds 99.1% of its correlation and the LFM
54.8%, which is the invariance that justifies the choice.

Three separate bugs came out of that one HFM measurement, each producing
a believable answer: render_real_doppler takes v/c and not 1 + v/c, so
the baseline was a 2x compressed pulse and ratios exceeded 100% -- which
a matched filter cannot do, and which is what exposed it; raw peaks
compared instead of correlation coefficients, so energy differences
masqueraded as mismatch; and a lag window starting at zero, which missed
the NEGATIVE peak shift an upsweep produces (v0.2's dt = -delta f_end/mu).

CRC-32 against its published check value 0xCBF43926, and RS(15,11) over
GF(16) correcting 2 symbol errors with 0 failures in 4000 trials. Beyond
t=2 the layers were measured TOGETHER: RS miscorrects 992 of 3000
four-error words to a valid-but-wrong codeword -- the minimum distance
d=5 showing, not a defect -- and the CRC catches all 992, none escaping.

Deferred and said so: carrier and chip-timing recovery, and an equaliser.
The demodulator is GIVEN the carrier phase. A few lines of arctangent
would not be a recovery loop.

181 cases, 95806 checks. Clean under gcc and clang in both precisions,
under ASan+UBSan, zero allocation by source audit and nm. Bellhop
cross-validation and generated-data reproducibility still pass.

## 0.12.0 — the turn rate measured, and the IMM in the tracker

A BRACKET DOES NOT DEGRADE GRACEFULLY.

v0.11's IMM reports a probability-weighted blend of its model turn rates,
so the estimate cannot leave [-omega, +omega]. With models set for
3 deg/s, against a five-state filter that estimates the rate instead:

    truth   IMM (bracketed)   CTRV (estimated)
    2 deg/s      1.27              1.02
    5 deg/s      2.62              4.99
    8 deg/s      1.82              5.64

The 8 deg/s row is the one that matters: the IMM reports LESS turn than
at 5. Once the truth leaves the bracket the models fit so badly that
probability drifts back to the constant-velocity model and the blend
collapses toward zero. That is what justifies carrying a fifth state,
and it is not visible from the 5 deg/s row alone.

The price is a nonlinear transition -- omega multiplies the velocity
terms -- so the covariance goes through a Jacobian rather than a constant
matrix. Both the transition and the Jacobian contain the same
sin(wT)/w and (1-cos wT)/w, and they switch to their series expansions at
the SAME threshold: a state and a covariance switching at different
points describe different filters. All six quantities are computed in one
place for that reason. The measurement Jacobian's fifth column is zero;
omega is observed only through where the target is predicted to be, which
is why several scans of turning are needed before it is worth anything.

AND THE TRAP IN DOING SO. The random walk on omega is the one parameter
with no counterpart in a constant-velocity filter, and setting it to zero
fails invisibly:

    q_w      reported sigma    estimate (truth 5.00 deg/s)
    0             0.11 deg/s         2.08
    0.005         0.91               4.97
    0.01          1.54               4.90
    0.02          2.67               4.89

With no process noise the covariance shrinks monotonically, the gain on
the fifth state dies, and the filter reports the SMALLEST uncertainty of
any setting about the MOST wrong answer. Confident and wrong is a worse
failure than uncertain and wrong, because nothing downstream can tell.
The default was 0.02 when the filter was written and is lowered to 0.01
on the strength of this sweep. That the reported sigma settles rather
than collapsing is correct, not a defect: omega can change at any moment,
and a filter that stops allowing for that has stopped being a tracker.

The fifth state costs 1.02x on a straight target over 12 runs -- cheaper
than expected, because the measurement never moves it.

Neither filter is better. CTRV measures a manoeuvre; an IMM reacts to
one, faster, because a model that already fits is waiting to take over.
Both ship, and the header says which to reach for.

imm_tracker_step() completes the integration v0.11 deliberately left out:
the same global-cost-ordered association and M-of-N management as
tracker_step(), gating on the COMBINED estimate rather than per model. A
per-model gate would let the worst-fitting model veto a measurement the
mixture accepts happily, and during a manoeuvre that is every model but
one; the mixture's covariance already widens when the models disagree,
which is exactly when the gate should be generous. Measured: two targets
turning opposite ways hold two established tracks through 40 scans with
the correct manoeuvre signs, and 367 false alarms over 200 scans produce
zero tracks, so M-of-N survives the swap to a mixture.

Also: wrap_pi moved from tracker.cpp's anonymous namespace to types.hpp
rather than being copied a third time.

172 cases, 43287 checks. Clean under gcc and clang in both precisions,
under ASan+UBSan, zero allocation by source audit and nm. Bellhop
cross-validation and the generated-data reproducibility checks still pass.

## 0.11.0 — real data, and a coefficient that was wrong for eleven releases

A PUBLISHED CHECK VALUE FOUND A BUG THAT ELEVEN RELEASES OF TESTING
DID NOT.

Until now the sound-speed equations were verified by mutual agreement:
Medwin, Mackenzie and Chen-Millero must agree inside their common
validity box, so a mistyped coefficient breaks the agreement. They agree
to about 0.1 m/s, so the check could never see an error worth 0.016 m/s
-- and there was one. chen_millero() held Wong & Zhu's (1995) ITS-90
coefficients for 40 of its 42 terms and Chen & Millero's 1977 originals
for A02 and A03. It was neither equation, and its header credited Chen &
Millero (1977) for something that was mostly Wong & Zhu.

Mutual agreement bounds how wrong you can be by how much your methods
differ. That is not the same as being right.

The fix was not to pick a number but to implement both equations. Against
the primary source, UNESCO Technical Papers in Marine Science 44 p. 48
and p. 50: check value 1731.9954 against a published 1731.995, and all
220 published table values within 0.0499 m/s -- the table is printed to
0.1, so 0.05 is exactly its rounding half-width. A test also proves the
two versions are legitimately different equations rather than one
equation and one typo: converting the temperature scale (t68 = 1.00024
t90) improves their agreement from 0.0208 to 0.0056 m/s, inside Wong &
Zhu's own stated revision size of 0.024.

The table is transcribed by hand from a scanned 1983 document with no
text layer, so tools/make_unesco_table.py verifies the transcription
against separately-entered coefficients and regenerates the header; CI
diffs the result, so a hand-edited value cannot survive.

Del Grosso (1974) added as a genuinely independent fourth equation --
different laboratory, different data, different form. It disagrees with
UNESCO by 0.41 m/s in the top kilometre and 3.93 m/s over the full
nominal validity box, the difference being that the box contains 26 C
water under 1000 bar, which no sea does. A validity box is a rectangle
in (S,T,P) and the ocean is not a rectangle. Worth reading against Chen
& Millero's own quoted standard deviation of 0.19 m/s: two independent
equations differ by twice the uncertainty either claims, which is the
real accuracy of "the speed of sound in seawater".

Six real WOA23 profiles, fetched over OPeNDAP one grid cell at a time so
every number's provenance is a URL that can be re-fetched. The Black Sea
earns its place: an 18.2 PSU surface makes the usual 35 PSU assumption a
20.3 m/s error, against at most 4.7 m/s anywhere else -- the same size as
the entire 24 m/s channel excess there, so the channel you trace is not
the one that exists. The shipped descriptions are tested against the
shipped data; an earlier draft called the Norwegian Sea upward-refracting
and put the North Atlantic axis at 1100 m, and neither survived contact
with the profiles.

A test that was measuring itself. Snell's invariant over real profiles
reports EXACTLY zero drift, which is suspicious rather than reassuring:
the tracer carries xi and reconstructs theta from it, so cos(theta)/c
largely measures whether acos and cos round-trip. Replaced with a
turning-depth check independent of how the state is stored -- at theta=0
the local sound speed must equal c_src/cos(theta_0). 186 turning points
across six real profiles, all within 1e-9 m/s. The first version of that
check reported errors up to 4.3 m/s that were entirely its own, because
the turn sits at exactly zero and zero is not > 0, so the point before it
also looked like a sign change.

An IMM filter, and an honest account of what it buys. Three models on
one four-state vector: constant velocity, and coordinated turns at
+/-omega, which are linear at known turn rate so no augmentation is
needed. Below |wT| < 1e-4 the transition switches to a series, since
sin(wT)/w is 0/0 at the limit; a test walks w down through the crossover
and gets 6.946e-N with no jump.

Against a single-model EKF tuned with process noise BETWEEN the IMM's
two, on a 3 deg/s turn:

    during the turn   37.31 m vs 36.68 m   0.98x
    after it ends     25.88 m vs 44.53 m   1.72x
    overall           44.91 m vs 50.74 m   1.13x

The gain is in RECOVERY, not in the turn -- during the manoeuvre the IMM
is marginally worse. That is not the usual description of an IMM and it
is what the measurement says. It costs nothing on a straight target
(27.40 vs 27.88 m), so the win is not more process noise in disguise,
and 1 - mu_CV is a manoeuvre detector that comes free, peaking at 0.991.
The IMM is NOT yet wired into tracker_step(); that is v0.12, and the
header says so.

166 cases, 43259 checks. Clean under gcc and clang in both precisions,
under ASan+UBSan, zero allocation by source audit and nm. Bellhop
cross-validation still passes.

## 0.10.0 — fusing what was already measured

Four quantities were produced by one subsystem and discarded at the
boundary of the next. No new physics; just joining them up.

Range rate into the tracker. The Doppler bank has estimated closing rate
since v0.3 and the filter ignored it, inferring velocity from position
history alone. The measurement is now 3-dimensional, with the chi-square
gate following the dimension -- a 3-dof quantile has no closed form, so
it is bisected on the exact CDF and reproduces 7.815 / 11.345.
Consistency holds: mean NIS 2.998 against a theoretical 3.

Radial-velocity error improves 3.45x at 3 scans, 3.19x at 6, and 0.99x
at 20. That last figure is published rather than cropped: by twenty
scans the filter's own estimate already beats the 2 m/s measurement, and
a measurement helps exactly as long as it beats the estimate you have.
The first version of the test asserted improvement at every step count,
failed, and was wrong about the physics rather than the code.

An unresolved Doppler bin reports 0 m/s, which is not a measurement that
the target is stationary. Fed in as one it reads a true 8.00 m/s closure
as 5.686. Hence a separate has_range_rate flag, and a test for it.

Ghost recognition, delivering what v0.8 wrongly expected of tracking.
v0.9 measured that time consistency cannot suppress a cross-template
ghost, because a ghost is exactly as consistent as the target. What
separates them is origin, not kinematics: the same arrival through a
different matched filter, so shared bearing, fixed offset, weaker, and a
DIFFERENT waveform label. The label check is the entire safety argument
-- two real targets lit by one sonar return the same waveform. A
line-astern formation is kept; the same geometry with different
waveforms is suppressed; the two runs differ only in that field.

Global-cost-ordered association. Every gating pair is built, sorted by
NIS, assigned best-first. Crossing targets keep their identities where
greedy-by-arrival-order swapped them. Still greedy, just over a global
ordering rather than arrival order, and not claimed to be optimal. The
gate was not touched: a swap is not a gating failure, since near the
crossing both measurements gate against both tracks perfectly well.

Forward-backward spatial smoothing, so MVDR survives the coherent
multipath v0.4 produces. Two coherent arrivals at +/-6 deg go from 4
spurious peaks to 2 correct ones. The cost is aperture -- resolution
falls from 7.16 to 11.46 deg -- and a test asserts that it does. The
first implementation had transposed indices and updated in place over
entries it had yet to read; it produced a plausible-looking covariance
that read a 9 deg source at 3.3. The invariant that catches it is
persymmetry.

Also: a false verification, and a guard for it. The four-configuration
sweep was run with -DPHANTOM_USE_FLOAT=ON. No such option exists; CMake
ignores unknown -D values silently, so both "float" builds were double
and all 41044 checks passed. A clean report for a run that never used
float is worse than a failure. CMakeLists.txt now rejects that name and
eight other near misses with a fatal error, the genuine float runs were
made, and the two builds are confirmed distinct by the Snell drift
(7.08e-08 against 1.67e-16).

count_established() added: Confirmed OR Coasting. A track that missed
the most recent scan is still a contact, and counting only Confirmed
reports the detector's instantaneous miss rate.

152 cases, 41044 checks. Clean under gcc and clang in both precisions,
under ASan+UBSan, zero allocation by source audit and nm. Bellhop
cross-validation still passes.

## 0.9.0 — tracking, and a roadmap claim corrected

A Pulse Descriptor Word has carried time, type and Doppler since v0.2 and a
bearing since v0.8. Nothing connected them across blocks, which is what turns a
list of detections into a picture.

The filter
  An EKF over a constant-velocity target: Cartesian dynamics because they are
  linear there, polar measurements because that is what a sonar produces, and
  all the nonlinearity confined to the measurement Jacobian. Discrete
  white-noise acceleration for Q, with sigma_a the single knob that sets how
  much manoeuvre the filter expects.

The check that matters
  A Kalman filter whose covariance is wrong still tracks. It just lies about
  how well, and then gates correct measurements out or accepts clutter, with
  nothing in its output to say so. The normalised innovation squared is
  chi-square with 2 degrees of freedom when the filter is consistent, so:

    mean NIS over 1980 samples   1.951   (expected 2.000)
    under the 95% gate           95.0%   (expected 95%)
    under the 99% gate           99.0%   (expected 99%)

  Every position test in the world passes on a filter that fails this one.

  The chi-square CDF with 2 dof is 1 - exp(-x/2), so gates invert in closed
  form and need no table: 5.991 at 95%, 9.210 at 99%.

What it is worth
  RMS position error 50.97 m raw against 20.89 m filtered, 2.44x better. And
  velocity, which a single detection cannot know at all: (4.17, -7.18) against
  a truth of (4.00, -7.00). The closing rate comes out at 7.683 m/s against
  7.486 -- the same quantity the Doppler bank measures directly, so the two
  subsystems are independent routes to one number.

Association and management
  Chi-square gating, greedy nearest-neighbour association, M-of-N confirmation,
  coasting through a miss and deletion after enough of them. Greedy is not
  optimal -- a global assignment does better when targets cross -- but it is
  O(tracks * measurements) with no allocation and its failure mode is a track
  swap, which is well understood rather than surprising.

  494 false alarms scattered over 300 scans produce ZERO confirmed tracks.

A roadmap claim, corrected in public
  The v0.8 roadmap said time consistency would finally suppress the
  cross-template ghosts that have been a documented limitation since v0.2.

  That claim was wrong. A ghost appears whenever the real arrival does, at a
  fixed offset set by the template cross-correlation, so it is EXACTLY as
  consistent over time as the target, moves with it, and forms its own
  perfectly healthy confirmed track. Nothing about its kinematics is
  objectionable.

  Measured: one target plus its ghost gives two confirmed tracks. The test
  asserts that outcome, with the reasoning in the test body, so the correction
  cannot quietly rot back into the optimistic version. Tracking kills false
  alarms -- which do not repeat. Ghosts repeat.

  Suppressing them needs the fixed offset and amplitude ratio recognised as a
  template artefact, which is a different mechanism and is now v0.10.

Two details that would have bitten later
  The bearing innovation is wrapped to [-pi, pi]. Without it a track near +/-180
  degrees produces a 2 pi innovation and diverges on its first update.

  (I - KH)P is symmetrised after every update. It is symmetric in exact
  arithmetic and drifts out of it in floating point; an asymmetric covariance
  eventually goes indefinite and the filter diverges with no warning at all.

41008 checks across 143 cases. Clean under GCC 15 and Clang 21 with -Werror,
under ASan+UBSan, and in both double and float builds. Zero allocation still
proven by nm. Bellhop cross-validation still passes.

## 0.8.0 — wideband beamforming, shading, MVDR

v0.7 built an array and then measured why it could not be used: phase steering
holds to 91 Hz of bandwidth off broadside, and this library transmits 12 kHz
chirps. This release fixes that and joins the array to the pulse analyser, so a
Pulse Descriptor Word finally carries a bearing.

Wideband, by steering per frequency bin
  A phase shift that is wrong across a band is exactly right WITHIN one bin. So
  transform, steer bin by bin, transform back -- an exact fractional delay
  rather than an interpolation of one.

  Two details decide whether it works. Bins above M/2 carry NEGATIVE
  frequencies, and steering them with the positive value breaks the spectrum's
  Hermitian symmetry so the inverse transform returns noise that still looks
  like a signal. And the forward model's sign must match the steering
  convention.

  Measured: a 5 ms chirp over 8-20 kHz at 35 degrees, 37x beyond the
  phase-steering limit, gives 17.0 dB of contrast between the on-target beam
  and the best off-target one -- and the matched filter peak lands at lag 400
  for a pulse placed at 400. A phase-steered beam would have smeared the chirp
  and moved it.

The structural check that earned its place immediately
  A narrowband tone pushed through the WIDEBAND path must reproduce
  array_factor in closed form. Worst departure across 21 steer angles: 5e-3.

  The first version of the test harness delayed the element at +x instead of
  advancing it -- the mirror of synthesize_plane_wave's convention, which v0.7
  had already validated against the Cramer-Rao bound. Every beam steered to the
  wrong side. The energy-contrast test only reported that the contrast was
  poor; the closed-form comparison said the pattern was mirrored, which is a
  different and far more useful sentence.

MVDR resolves what the aperture cannot
  Conventional beamforming cannot separate two sources closer than one
  null-to-peak spacing WHATEVER the SNR: resolution is set by the aperture and
  nothing else. MVDR places nulls instead of scanning a fixed pattern.

  On a 16-element half-wave array whose conventional limit is 7.16 degrees, two
  incoherent sources 4.30 degrees apart give one conventional peak and two MVDR
  peaks, at +/-2.145 against a truth of +/-2.149. And MVDR agrees with
  conventional on a single source (11.000 vs 10.988 for an 11 degree truth), so
  the resolution is not bought with a bias.

  Solved by complex Cholesky rather than an explicit inverse: R is Hermitian
  positive definite after loading, so one forward/back substitution per angle.

  Diagonal loading is verified as NECESSARY, not decorative. A rank-1
  covariance over 16 elements fails the factorisation outright and mvdr_power
  returns 0. Reporting the failure beats returning noise shaped like an answer.

Shading
  The shaded pattern is the DTFT of the window, so its sidelobe level IS the
  window's. Measured against the published values: -13.26 / -31.5 / -42.7 /
  -58.1 dB for uniform / Hann / Hamming / Blackman, with mainlobe widths of
  1.00 / 1.64 / 1.47 / 1.90x and gain given up of 0 / -1.76 / -1.34 / -2.37 dB.
  Every window widens the beam and costs gain -- asserted, not just observed.

The integration
  The analyser run per beam: a source at -22 degrees peaks in the -20 degree
  beam of a 5 degree scan, at 38 dB. A PDW now carries time, type, Doppler AND
  bearing.

  The two-stage API exists because the memory-versus-time trade is real and
  belongs to the caller: transforming every element once costs N*M complex
  values, re-transforming per beam costs 1464 FFTs instead of 24 for a
  24-element array and 61 beams.

Scope kept honest, in docs/validation.md 15
  MVDR needs incoherent sources -- the coherent multipath this same library
  produces in v0.4 defeats it, and spatial smoothing is not implemented. No
  tracking: a PDW carries a bearing but nothing associates detections across
  blocks. Uniform line arrays only, no element directivity or calibration
  errors. Reverberation is still a scalar-beamwidth sonar equation, not a
  per-beam field.

40964 checks across 131 cases. Clean under GCC 15 and Clang 21 with -Werror,
under ASan+UBSan, and in both double and float builds. Zero allocation still
proven by nm. Bellhop cross-validation still passes.

## 0.7.0 — line arrays, beamforming and bearing

The last purely structural gap. Everything before this release was
single-channel: a hydrophone knows WHEN a ping arrived and what shape it was,
but not where from -- and v0.6 showed reverberation scales with the ensonified
area, of which the azimuthal half can only be shrunk with a beam.

The bearing bound, and why it is the spatial twin of v0.2
  Same derivation as the arrival-time bound, with the element index in place of
  time and Sum n'^2 = N(N^2-1)/12 about the array centre in place of the
  waveform's mean-square bandwidth:

    var(theta) >= 6 / (rho (k d cos theta)^2 N (N^2 - 1))

  Two properties fall straight out and both are verified. Accuracy improves as
  N^(-3/2) rather than N^(-1/2), because elements buy signal AND aperture --
  doubling N should gain 2^1.5 = 2.828, measured 2.845/2.833/2.830 across
  N = 8 to 64. And it degrades as 1/cos(theta), reproduced to 1e-6 out to 75
  degrees, because a target at endfire sees no projected aperture.

  A 32-element half-wave array is bounded at 0.247 degrees against a 3.17
  degree beamwidth -- one thirteenth of a beam. Bearing comes from phase across
  the aperture, not from which beam lit up.

  The estimator lands at 1.13/1.02/0.99 x the bound across 14 dB of SNR. No
  factor-of-four surprise this time: unlike v0.2's arrival-time case, where the
  coherent and envelope bounds differ by 4.16 and picking the wrong one made an
  efficient estimator look broken, here there is no carrier to lose. The phase
  gradient across the aperture IS the signal.

Beam pattern, pinned by closed forms
  A pattern that is subtly wrong still looks like a beam, so: unity at the
  steer angle with nothing exceeding it over +/-90 degrees; nulls exactly at
  sin(theta) - sin(theta_0) = m lambda/(N d); and a first sidelobe converging
  to the -13.26 dB signature of uniform shading -- -12.797 at N=8, -13.260 at
  N=128.

  Steered beams get the exact null position, asin(sin(theta_0) + lambda/Nd)
  minus theta_0, not the small-angle form. A beam at 60 degrees is genuinely
  broader because the projected aperture shrinks as cos.

Array gain, measured rather than asserted
  10 log10(N) through the actual beamformer over 400 trials: 5.92/12.07/18.36
  dB for N = 4/16/64 against 6.02/12.04/18.06 predicted. The formula is
  trivial; the beamformer is not, which is why the Monte Carlo is there.

Split-beam
  Exact to 1e-16 across the mainlobe, unambiguous precisely out to the first
  null because |dphi| = pi there, and demonstrably wrapping beyond it -- a
  source at 8.96 degrees reads as 1.76. Documented limit, not a defect.

It joins up with v0.6
  Ten times the elements is a tenth of the beamwidth is exactly 10.00 dB of
  echo-to-reverberation ratio -- the same 10 dB per decade of ensonified area
  that v0.6 charges. The same change also buys 10.0 dB of array gain against
  isotropic noise. Two independent mechanisms agreeing is a stronger check than
  either alone.

A limit stated rather than hidden
  Phase steering only holds while the signal stays correlated across the
  array's traversal time: 249 Hz of usable bandwidth at 15 degrees, 91 Hz at
  45. This library's own waveforms are 12 kHz chirps, so a phase-steered array
  cannot handle them off broadside. narrowband_bandwidth_limit_hz returns the
  number instead of the beamformer quietly smearing the beam.

Also fixed
  bearing_crlb_rad guarded endfire with cos > 1e-9. cos(pi/2) is 6e-17 in
  double but -4e-8 in float, so the guard stopped firing in a single-precision
  build and the bound came back finite where there is no bearing information at
  all. Same failure and same fix as spreading_loss_db in v0.4; the comment now
  cites that precedent.

Scope kept honest, in docs/validation.md 14
  Narrowband only. Uniform line arrays only -- no shading, no planar geometry,
  no element directivity. One source: the estimator finds a peak, and two
  sources inside a beamwidth are not resolved. And the array is NOT wired into
  the analyser yet -- beamforming and pulse analysis are still separate pieces,
  which is v0.8.

40915 checks across 123 cases. Clean under GCC 15 and Clang 21 with -Werror,
under ASan+UBSan, and in both double and float builds. Zero allocation still
proven by nm. Bellhop cross-validation still passes.

## 0.6.0 — reverberation

v0.5 capped the surface scattering loss and said the energy "goes into a
diffuse field the ray model does not carry". This is that field coming back.
At short range it, not the ambient noise the analyser assumes, is what an
active sonar competes against -- and that changes which knobs help.

The result this module exists for
  Write the echo and the reverberation with the same source level and the same
  two-way loss:

    EL      = SL - 2 TL + TS
    RL      = SL - 2 TL + S_s + 10 log10(A)
    EL - RL = TS - S_s - 10 log10(A)

  SL and TL cancel EXACTLY. Verified over nine combinations of source level
  (160/200/240 dB) and transmission loss (40/66/90 dB): the ratio comes out
  5.2288 dB every time, identical to four decimals. Against reverberation a
  louder transmitter raises the target and the background together and buys
  nothing at all.

  What does help is shrinking the ensonified area: a tenth of the pulse length
  or a tenth of the beamwidth, each worth exactly +10.00 dB. And a chirp's
  effective pulse length is 1/B rather than its duration, so compressing 20 ms
  to its 83 us resolution is worth +23.8 dB -- the time-bandwidth product in dB.
  That is the design argument for pulse compression, and it is why the
  analyser's waveforms have been chirps since v0.2.

Two decay laws that identify the mechanism
  Boundary reverberation falls as 30 log10(r): two-way spreading costs 40, the
  growing annulus gives back 10. Volume reverberation falls as 20, because the
  shell grows as r^2. Both reproduced exactly across three decades. A measured
  decay of 20 rather than 30 says the water column is scattering, not the
  boundary.

Scattering strengths
  Lambert (mu + 20 log10 sin theta) for the bottom, reproducing its closed form
  exactly and returning mu at normal incidence. Chapman-Harris for the
  wind-driven surface.

  Chapman-Harris is flagged in docs/validation.md 11 rather than trusted: the
  formula is written out in the header so a reader can check it, and its
  BEHAVIOUR is verified -- monotone in wind over 3-40 knots, monotone in angle
  over 2-80 degrees, correct structure at exactly 30 degrees where the angular
  term vanishes -- but the coefficients have not been checked against the 1962
  paper. Two reference values in this repository have already turned out to be
  misremembered rather than wrong in the code (Thorp at 50 kHz in v0.4, a
  critical angle in v0.5), so this one is marked rather than assumed.

Reverberation-limited range
  Reverberation falls as 30 log10(r); ambient noise does not fall at all. Their
  crossover is where a bigger transmitter starts helping again: 24.7 km at a
  40 dB ambient, 5.3 km at 60, 1.1 km at 80. Verified to actually be a
  crossover -- the two levels agree to 0.1 dB at the returned range.

What it does to the detector
  Reverberation decays tens of dB across a single processing block, so no fixed
  threshold is right anywhere. Measured on a block whose background falls
  21.6 dB end to end: CA-CFAR gets 20/20 detections of a target buried in the
  decayed region with 3 false alarms across 20 empty blocks, while a single
  threshold at the block mean sits 9.5 dB below the near field and 12.2 dB
  above the far field. That gap is precisely what CFAR closes, and until now
  the test suite had never made it visible.

Scope kept honest
  Reverberation is a LEVEL, not a field: the envelope generator gives the right
  level and the right post-correlation statistics, but a true series is the
  scatterer field convolved with the transmitted pulse and so is correlated
  over the pulse length, where white noise scaled by an envelope is correlated
  over one sample. Documented at the function and in validation.md 13.

  Also: reverberation comes from the sonar equation with spherical spreading,
  not from the traced eigenrays; no bistatic case; Lambert carries no frequency
  dependence.

39926 checks across 113 cases. Clean under GCC 15 and Clang 21 with -Werror,
under ASan+UBSan, and in both double and float builds. Zero allocation still
proven by nm. Bellhop cross-validation still passes.

## 0.5.0 — boundary reflection losses

Until now every surface and bottom reflection was perfect. In shallow water
most paths bounce, so this was the largest remaining overstatement in the
transmission loss: a four-bounce path was reported as a peer of the direct one.

Bottom — the Rayleigh coefficient
  A fluid-fluid interface, pinned by three independent closed-form limits:

    critical angle      arccos(c1/c2), from the speed ratio alone
    normal incidence    (Z2-Z1)/(Z2+Z1), the impedance ratio
    below critical      |R| = 1 exactly, lossless total internal reflection

  Any one of those could be reproduced by a wrong expression; all three
  together could not. 24.62 degrees for medium sand, and that single number
  decides which paths survive to long range in shallow water. A bottom SLOWER
  than the water has no critical angle at all and leaks everywhere, which is
  why mud is acoustically so much worse than sand.

Sediment attenuation, and why it is not optional
  Published as dB per wavelength; enters as a complex sediment speed via
  delta = alpha / 54.575. |R| < 1 then holds at EVERY angle including below
  critical. Without it, sub-critical rays are trapped forever and shallow-water
  range is unbounded, which is not what the ocean does. The default sand costs
  0.81 dB per bounce at half the critical angle -- 8 dB over ten.

  The complex square root has two branches and only one gives a wave that
  decays into the sediment. Rather than reason about sign conventions, the
  implementation takes the branch that conserves energy: the other yields
  |R| > 1, which no passive interface can do. A sweep over 5 speeds x 4
  densities x 5 attenuations x 91 angles -- 9100 combinations -- gives
  max |R| = 1.000000000000.

Surface — and a cap that is the honest part
  Coherent scattering from a rough pressure-release surface,
  |R| = exp(-G^2/2) with G = 2 k sigma sin(theta), and Pierson-Moskowitz
  H_1/3 = 0.0246 U^2 for wind.

  At 0.5 m seas, 10 kHz and 20 degrees grazing this gives a 700 dB loss. That
  is arithmetic, not physics: the scattered energy is not destroyed, it goes
  into a diffuse field a ray model does not carry. The path is not 700 dB down,
  its SPECULAR part is. So the reported loss is capped -- 30 dB by default --
  and the header says why, rather than deleting a path that is still in the
  water. Filling that gap properly is reverberation, and it is v0.6.

Applied along the path
  The grazing angle at each bounce comes from the Snell invariant rather than
  from having tracked it: cos(theta) = xi * c at the boundary depth, exact in a
  range-independent ocean. Verified against a hand computation -- 2 surface + 3
  bottom bounces at 10 degrees gives 6.7532 dB, matching the sum of the
  individual losses.

  The effect, for a 200 m duct at 3 km, 5 kHz, 8 m/s wind, sand bottom:

    srf btm  graze   TL before   TL now
      0   0      -       71.80    71.80
      1   1   6.83d      70.67    82.85
      2   2  13.93d      70.95   132.73
      4   3  24.69d      71.62   199.51

  Paths that all sat within 1 dB of each other now span 128 dB, and the direct
  path is unchanged as it must be.

Also fixed
  A hardcoded "expected" critical angle of 14.62 degrees for 1550 m/s sediment
  was wrong -- the correct value is 14.593, which the implementation had. Same
  failure mode as the Thorp reference in v0.4: my arithmetic, not the code's.

Scope kept honest, in docs/validation.md 12
  No diffuse field. Plane-wave flat-interface reflection, so no beam
  displacement near the critical angle, no sediment layering, no shear.
  Caustics still flagged rather than computed. Monostatic multipath.
  Range-independent.

15642 checks across 104 cases (the count jumps because of the 9100-combination
energy sweep). Clean under GCC 15 and Clang 21 with -Werror, under ASan+UBSan,
and in both double and float builds. Zero allocation still proven by nm.
Bellhop cross-validation still passes.

## 0.4.0 — ray-acoustics coupling: eigenrays, transmission loss, multipath

The two halves of the library used to share only a sound speed. Now the ray
tracer supplies the delays and levels the echo synthesiser transmits, so a
multipath arrival structure comes out of the physics rather than out of a
parameter.

Eigenray search, exact rather than interpolated
  The tracer already honours a range budget to machine precision, so setting
  that budget to the receiver range makes the ray's final point THE arrival --
  no polyline resampling and none of the chord-versus-arc error that cost 12 m
  of apparent disagreement in the Bellhop work. What remains is a root find on
  z(theta_0) - z_receiver.

  Isovelocity check: launch angle, path length, travel time and Jacobian all
  reproduce the straight-line geometry (6.842773 deg vs 6.842773, 5035.871 m vs
  5035.871). A 200 m duct at 3 km yields 14 eigenrays with bounce counts from
  0/0 to 4/3, ordered and distinct.

Exact arc length in the tracer
  Absorption is quoted per unit path, not per unit range, so TraceResult now
  carries path_length_m accumulated exactly: R|dtheta| on each constant-gradient
  segment, |dz/sin| on each isovelocity one. No approximation.

Geometric spreading from the ray tube
  Conserving power in the tube between neighbouring rays:

    TL = -10 log10[ c_rcv cos(theta_0) / (c_src r cos(theta_rcv) |dz/dtheta_0|) ]

  The cos(theta_rcv) is the tube's cross-section perpendicular to the ray;
  leaving it out costs 10 log10(cos theta), which is how I found it -- the first
  derivation was 3 dB off spherical spreading at zero angle.

  In isovelocity water this must reduce to 20 log10(R). It does, to 0.0000 dB
  across nine geometries from 1 to 20 km, with the Jacobian obtained by finite
  difference THROUGH the tracer -- so the agreement exercises the search, the
  arc-length accumulator and the formula at once.

Thorp absorption
  0.069 dB/km at 1 kHz, 1.19 at 10 kHz, 34.1 at 100 kHz. Over 3 km at 100 kHz
  absorption exceeds spreading by 30 dB, which is the whole reason long-range
  sonar is low-frequency.

  My first test used 15.9 dB/km at 50 kHz as the "published" value and failed.
  An independent evaluation of the published coefficients gives 17.47, matching
  the implementation exactly -- the reference value was misremembered, not the
  code. Tolerances now sit at the ~10% the fit is actually worth.

Caustics flagged, not answered
  Where the tube collapses, ray theory predicts infinite intensity. That is a
  failure of the method, not a property of the ocean, so such paths are flagged
  and left without a level rather than given a number that looks like one.

A precision limit worth knowing about
  depth_tolerance_m DISCARDS roots it cannot polish, so a tolerance tighter than
  the build can achieve removes paths silently instead of degrading them.
  Measured on the 200 m duct, 14 paths present:

    tolerance   double  float
      0.001 m      14      0
      0.010 m      14      3
      0.100 m      14     14

  The default now follows the build. This is exactly the kind of failure that
  looks like "the ocean has fewer paths than I thought".

Fixed along the way
  - spreading_loss_db guarded a vertical arrival with cos > 0. cos(pi/2) is
    6e-17 in double and -4e-8 in float, never zero, so it returned a few hundred
    dB instead of declining to answer. Now an epsilon.
  - EchoSpec gained extra_delay_s as its SECOND member, and every positional
    aggregate initialiser silently became a different echo: a -8 dB target
    strength was read as a -8 SECOND delay and the ghost was dropped for
    arriving before the ping. Tests and examples now use designated
    initialisers, so the next insertion is a no-op instead of a puzzle.

Scope kept honest
  Ray theory with its caustics. No bottom loss -- reflection is specular and
  lossless, so a bottom-bounced path is reported louder than it is. Monostatic
  only. No surface scattering or bubble loss. Range-independent. Target strength
  is still a dB figure, not a scattering model. All listed in
  docs/validation.md 11 rather than left to be discovered.

6440 checks across 92 cases. Clean under GCC 15 and Clang 21 with -Werror,
under ASan+UBSan, and in both double and float builds. Zero allocation still
proven by nm. Bellhop cross-validation still passes.

## 0.3.0 — Doppler bank and echo synthesis

The countermeasure half, and the Doppler coverage the analyser needed to
support it. Both are verified by closed loop: a ping is detected, an echo is
synthesised from the resulting descriptor, and the echo is fed back through the
analyser. A shared sign error between the two halves then shows up as a round
trip that does not close, rather than as two tests that agree and are both
wrong.

Doppler bank
  Bin spacing is derived per waveform family and checked against the measured
  matched filter loss, not assumed:

    CW    output is sinc(delta f0 T); second-order expansion
    LFM   a mismatched chirp leaves a residual quadratic phase. The detector
          absorbs its least-squares linear part as the range bias of 6.2; what
          remains has RMS 2 pi TB delta / sqrt(180), so delta_1dB ~ 1/TB
    HFM   time scaling maps the family onto itself EXACTLY -- the phase
          difference is constant when v = (alpha-1)/alpha -- so only the
          duration mismatch remains, an amplitude factor of (1 - |delta|)

  Requested 1 dB gives 1.08 dB measured for the LFM, 1.55 dB for the HFM. Over
  +/- 20 m/s at 1 dB, the same 8-20 kHz 20 ms sweep needs 14 bins as a CW, 5 as
  an LFM and 2 as an HFM. Each bin is another correlation per block and another
  128 kB of spectrum, so that ratio is the design argument for hyperbolic
  sweeps stated as a cost.

  Documented as a small-mismatch expansion. At 15 m/s on that LFM it predicts
  49 dB where 9.9 dB is measured, because the true loss saturates while the
  expression runs away. It spaces bins; it does not model deep mismatch.

Echo synthesis
  Delay and Doppler are both TWO-way: dt = 2dr/c and alpha = (c+v)/(c-v), the
  exact form rather than 1 + 2v/c (0.4% apart at 30 m/s, four samples across a
  20 ms pulse). The factor of two in each is the likeliest bug in an echo path,
  so the round trip asserts the reading is closer to the two-way scale than to
  the ghost's own velocity. Ghost swarms, target strength and extended targets
  all round-trip to within a sample and 6% of amplitude.

Anti-phase cancellation: quantified, not implemented
  20 dB of cancellation needs the timing held to 1.33 us -- 1.99 mm of path at
  1500 m/s. Bandwidth is nearly irrelevant: with a 2 us error at 12 kHz,
  widening the band from 0 to 12 kHz moves the residual by 0.34 dB. Past a
  quarter period the canceller ADDS up to 6 dB.

  The header's first draft claimed a bandwidth limit. Measuring it showed the
  limit is timing, and the claim was corrected rather than kept because it
  sounded right. Two further effects are named and left unmodelled, both of
  which make the real figure worse: distributed hull scattering, and the
  latency floor on tau. Generating false targets is the achievable
  countermeasure; cancelling the real one is not.

The two halves cross-check each other
  In countermeasure_loop an 8 m/s ping lands in the 12 m/s bin, and the 4 m/s
  residual produces an arrival-time error of +0.0885 ms. The wideband coupling
  formula from v0.2 predicts dt = -(v/c) f_end/mu = +0.0885 ms. Two results
  derived independently, agreeing to four decimals -- and CI now asserts both
  numbers, so they cannot drift apart silently.

Fixed along the way
  - add_doppler_bank rounded the interval count DOWN, so the realised spacing
    could be twice the design spacing and the straddling-loss guarantee did not
    hold. Now rounds up.
  - The detector reported whichever cell the dead-time window happened to end
    on when the window's maximum sat on its trailing edge: a one-sample bias in
    the good case and a wholly wrong arrival in the bad one. It now climbs to
    the local maximum.
  - PulseBank is megabytes and a 64-template one overflowed the stack in a
    test. Documented at the type, with the size formula.

Scope kept honest
  A bank only resolves what it spans: an echo stretched to 28 ms is detected by
  a 20 ms template and mis-reported in both time and type. countermeasure_loop
  transmits one deliberately and says so. Echo levels are specified in dB, not
  computed from a scattering model, and there is no transmission loss along the
  traced path -- the two halves still meet only through the sound speed. That
  coupling is v0.4.

6226 checks across 82 cases. Clean under GCC 15 and Clang 21 with -Werror,
under ASan+UBSan, and in both double and float builds. Zero allocation still
proven by nm. Bellhop cross-validation still passes.

## 0.2.0 — PingAnalyzer: matched filter bank, CFAR detection, PDW output

The stage the original specification skipped. An echo synthesiser with no
input is not a module: detecting and characterising the incoming ping is the
harder half of a countermeasure and the half with the latency budget, so it
ships before the synthesiser.

New
  fft.hpp/cpp             radix-2, in-place, caller-owned twiddles
  waveform.hpp/cpp        CW / LFM / HFM, closed-form phase, true Doppler
  matched_filter.hpp/cpp  FFT overlap-save replica correlation
  ping_analyzer.hpp/cpp   filter bank + CA-CFAR + Pulse Descriptor Words
  examples/ping_intercept streaming detection over 2 s, scored against truth

Zero dependencies still means zero: the transform is written here rather than
pulled in, and nm over the built archive shows the library referencing only
libm and memset.

Verification -- 6106 checks, every one against a closed form or a bound
  FFT vs a direct O(N^2) DFT sharing no code      2.3e-16 relative
  FFT correlation vs direct O(N*L) correlation    1.1e-15 relative
  d(phase)/dt vs the specified f(t), LFM and HFM  1.9e-11 relative
  matched filter peak / width / coherent gain     A*E/2, 1/B, L/2 -- all met
  classification, 4 waveforms at -4.4 dB SNR      100/100

Arrival time vs the Cramer-Rao bound
  0.93-1.02x the bound across 20 dB of SNR, flat -- the signature of an
  estimator limited by noise rather than by systematic interpolation error.

  Getting the comparison right mattered more than the estimator did. A
  magnitude detector is bounded by the RMS bandwidth about the CENTRE
  frequency; only a carrier-phase-coherent receiver is bounded by f_rms about
  zero. For an 8-20 kHz sweep those differ by 4.16x, and checking an envelope
  detector against the coherent bound made an efficient estimator look four
  times worse than it was. Both bounds are now computed and the applicable one
  named. The Fisher information is derived two independent ways, agreeing to
  0.06%.

A correction to the textbook LFM range-Doppler formula
  Underwater, Doppler is a time SCALING, not a frequency shift -- v/c is ~500x
  larger than for an airborne radar at the same speed. Carrying that through,
  the delay bias uses the sweep's END frequency, not its centre:

      dt = -(v/c) * f_end / mu

  Verified on upsweeps, downsweeps and narrow sweeps: the wideband form agrees
  within 3 samples where the narrowband one is off by 13 on the 8-16 kHz case.
  Derivation in docs/math_spec.md 6.2.

  Same physics measured across the three waveform families, 60 ms 8-20 kHz:
  at 30 m/s closing the LFM loses 13.1 dB and the HFM 0.30 dB. That is why
  HFM dominates underwater -- a time-scaled HFM is a delayed HFM.

A silent failure mode, reproduced and then fixed
  CA-CFAR estimates noise from cells either side of a guard band. If the guard
  is narrower than the target's response, the target leaks into its own
  training cells and raises its own threshold: the stronger the pulse, the
  higher the bar it must clear. Nothing is reported and the detector looks
  healthy.

  A chirp compresses to ~fs/B cells and survives a small guard. A CW does not
  compress at all. Measured with a 32-cell guard against a 1920-sample CW
  response: CW 0/20 detected, LFM 20/20. Three of four waveforms still working
  is exactly what makes this hide. suggested_cfar_guard() and
  suggested_dead_time_s() now return the right values, and the test pins both
  the failure and the fix.

Performance
  fft_forward()  8192 points   307 us   (a 3x win from hoisting the inverse
                                         flag to a template parameter: as a
                                         runtime bool it cost a branch per
                                         butterfly and blocked vectorisation)
  analyze_block()  4 templates 2.26 ms
  29x real time on one core at 96 kHz; 65 ms detection latency, set by the
  block length rather than the arithmetic.

  The transform is a plain radix-2 with no SIMD, roughly an order of magnitude
  off a vendor-tuned one. That is the price of the dependency promise, and it
  is isolated behind two functions if you would rather pay the dependency.

Documented limitations
  Cross-template ghosts: a bank whose waveforms share a band reports one
  arrival more than once -- an LFM arrival lights the HFM template 10.5 dB
  down at a shifted lag. Correct behaviour for a matched filter bank;
  suppressing it needs association logic across detections, which is a
  tracker's job and v0.4 work. Templates are zero-Doppler, so a fast target is
  detected rather than matched. No bearing: that needs an array.

Clean under GCC 15 and Clang 21 with -Werror, under ASan+UBSan, and in both
double and float builds. Bellhop cross-validation still passes.

## 0.1.1 — cross-validate the ray tracer against Bellhop

Closes the one gap v0.1 shipped with. Everything the test suite checked was
internal — closed forms and invariants — which catches implementation errors
but not a shared misunderstanding of the physics. This is the external check,
against Bellhop from the Acoustics Toolbox (at_2026_7, gfortran 15.3).

Design
  Both codes read the SAME .env file: phantom_trace parses Bellhop's own input
  format rather than a transcribed copy, because a transcription bug would
  look exactly like agreement. Interpolation is pinned to C-linear on both
  sides — phantom_trace refuses to run on a PCHIP or spline file, since that
  would measure the interpolation scheme rather than the tracer. Beam counts
  are explicit, so both codes launch the same rays.

Result — turning depth vs Bellhop's integration step
  The turning depth is the metric worth quoting: one number per ray, no
  interpolation anywhere, exact here by construction.

    step 500 m -> 0.084   m        step  50 m -> 0.0017  m
    step 200 m -> 0.014   m        step  20 m -> 0.00036 m
    step 100 m -> 0.0024  m        step  10 m -> 0.00014 m

  0.14 mm over a 101 km path. The direction matters as much as the magnitude:
  Bellhop converges toward this library's answer, which is what should happen
  if the analytic arc solution is the exact limit of a stepping integrator.

Full-path RMS, and why it is a ceiling rather than a measurement
  0.274 m for trapped rays, 0.635 m with 62 surface/bottom bounces, at
  Bellhop's default step. Those numbers are limited by THIS COMPARISON, not by
  either tracer: the residual halves exactly every time the library's output
  sampling doubles (1.055 m at x8 -> 0.0068 m at x512), which is first order,
  the signature of chording across the parabolic nose at each turning point.
  The tracer emits one point per layer crossing, so a 27-point profile yields
  a ray sampled every few kilometres.

  Reporting the RMS without that sweep would present a limitation of the
  measurement as a property of the code. It is published as Case D.

  --refine N splits each tabulated layer into N sublayers. This does not
  change the modelled ocean — c(z) is linear within a layer, so inserted
  points lie on the same line, and ray_is_invariant_under_layer_refinement
  already pins that the traced answer does not move.

Added
  tools/bellhop_compare/setup_bellhop.sh   fetch + build the toolbox. Two
      Linux fixes over the shipped Makefile: drop the macOS-only -Wa,-q, and
      build serially — undeclared Fortran .mod dependencies make `make -j`
      race and fail.
  tools/bellhop_compare/phantom_trace.cpp  .env parser + tracer -> CSV
  tools/bellhop_compare/make_env.py        cases from the canonical Munk table
  tools/bellhop_compare/compare.py         runs both, four cases, figure, JSON
  CI job with the toolbox build cached on the setup script's hash

Licensing: the Acoustics Toolbox is GPL-3. It is fetched onto the user's
machine and invoked as a separate process exchanging text files; nothing is
linked, vendored or redistributed, so no GPL obligation attaches here.

Scope: ray geometry only, range-independent profiles only. Says nothing about
transmission loss, beam amplitudes or caustics, none of which v0.1 computes.

## 0.1.0 — ocean model and analytic ray tracer

First release. Scope is deliberately one engine done properly rather than
three sketched: the ocean model and the ray tracer, verified against closed
forms, with the gaps documented rather than papered over.

Sound speed
  Three equations, because the usual one is mislabelled: the six-term
  polynomial widely circulated as "the UNESCO equation" is Medwin's
  simplification, valid only to 1000 m — at or above the deep sound channel
  axis, so it is the wrong tool for SOFAR work. Implements Medwin,
  Mackenzie (default, valid to 8000 m) and Chen & Millero (the actual UNESCO
  algorithm, pressure-based), plus Leroy-Parthiot depth->pressure and the
  Munk canonical profile. The three are cross-validated against each other,
  which catches any mistyped coefficient.

Ray tracer
  Analytic, not stepped. The profile is piecewise linear in depth, so within
  a layer the ray is exactly a circular arc: one closed-form step per layer
  crossing, no integration error, no step size to tune. Travel-time budgets
  invert the arc in closed form, so both range and time limits are honoured
  to machine precision instead of rounding back to a layer boundary.

Verification (791 checks, all against closed forms or invariants)
  Snell invariant over 29 rays x 100 km      1.7e-16 relative drift
  ray path vs analytic circle                1.9e-16 relative radial error
  layer refinement, 2 vs 1001 points, 40 km  4.6e-12 m
  range / time budgets                       1e-8 m / 1e-12 s
  Mackenzie vs Chen-Millero                  0.53 m/s max disagreement

  Clean under GCC 15 and Clang 21 with -Wconversion -Wsign-conversion
  -Wold-style-cast -Wdouble-promotion -Werror, and under ASan+UBSan. Passes
  in both double and float builds with precision-aware tolerances.

  Zero-allocation is proven, not asserted: nm over the built archive shows
  the library references only libm and memset.

Restated from the original specification
  The "< 500 ns per processing cycle" requirement was undefined — a 100 km
  trace crosses ~1100 layers and cannot fit that budget on any hardware. The
  budget is now stated per arc step, the unit that actually constrains a
  control loop, and measures 52 ns.

Documented limitations
  No external reference comparison yet — Bellhop cross-validation is the top
  v0.2 item, and until it runs, agreement is expected rather than
  demonstrated. Geometric acoustics only, range-independent, 2-D. Shadow
  fractions depend on ray density and are still drifting ~0.5 points per
  doubling at 5761 rays; the convergence table is published rather than a
  single flattering number.

Also: Apache-2.0 (patent grant), EXPORT_NOTICE.md, CI across three platforms
with source and link-level allocation audits, Cortex-M7 toolchain file, and
docs/hardware.md specifying the tank bench — including why it must run at
40-200 kHz rather than at realistic sonar frequencies.
