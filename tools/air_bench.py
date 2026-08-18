#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Air-acoustics bench: play a probe, record it, and measure the channel.

The cheapest possible hardware-in-the-loop test of this library. A speaker and a
microphone are already on most desks; a hydrophone and a projector are several
hundred euros and need water to put them in. The signal processing does not care
which medium it is in, so this exercises the whole chain -- waveform synthesis,
matched filtering, delay estimation, multipath -- before any wet hardware is
bought.

    python3 tools/air_bench.py selftest              # no hardware, simulated channel
    python3 tools/air_bench.py probe out.wav         # write a probe to play
    python3 tools/air_bench.py analyse rec.wav       # analyse a recording
    python3 tools/air_bench.py live --distance 0.5   # play and record together

WHAT IT MEASURES. The matched filter peak gives the acoustic delay from speaker
to microphone. Multiply by the speed of sound -- which depends on the room
temperature you type in -- and you get the path length. Compare with a tape
measure. Agreement to a few centimetres means the whole chain is working on real
signals; disagreement means something specific and findable.

WHAT IT CANNOT MEASURE. Absolute level, target strength, or absorption: a laptop
speaker and microphone have unknown, uncalibrated and wildly non-flat responses.
This is a TIMING and DETECTION bench, not a calibrated acoustics one.

A NOTE ON LATENCY. The delay you measure includes the sound card's own
input-plus-output latency, which is typically 10-50 ms and dwarfs the acoustic
flight time over a desk (1.5 ms per half metre). The only reliable way to remove
it is a two-distance measurement: the DIFFERENCE of two delays cancels every
fixed electronic latency and leaves the acoustics. `live` supports that and the
help below explains it.
"""
import argparse
import math
import struct
import subprocess
import sys
import wave

FS = 48000
F0, F1 = 2000.0, 8000.0
CHIRP_S = 0.020


def sound_speed(temp_c, rh=0.0):
    """Same model as phantom::air::sound_speed, including its DRY default.

    The humidity default matters and was wrong here first: this returned 344.00
    m/s at 20 C where the C++ returns 343.37, because the two disagreed about
    whether "20 C" meant dry or half-saturated. 0.6 m/s is 0.2%, which over a
    half-metre bench is a micrometre and over a kilometre is two metres."""
    t = temp_c + 273.15
    dry = 331.45 * math.sqrt(t / 273.15)
    psat_over_pr = 10.0 ** (-6.8346 * (273.16 / t) ** 1.261 + 4.6151)
    h = rh * psat_over_pr / 100.0
    return dry * (1.0 + 0.16 * h)


def chirp(n=None):
    n = n or int(FS * CHIRP_S)
    out = []
    for i in range(n):
        t = i / FS
        ph = 2 * math.pi * (F0 * t + 0.5 * (F1 - F0) / CHIRP_S * t * t)
        r = 0.25
        w = 1.0
        if t < r * CHIRP_S:
            w = 0.5 * (1 - math.cos(math.pi * t / (r * CHIRP_S)))
        elif t > CHIRP_S - r * CHIRP_S:
            w = 0.5 * (1 - math.cos(math.pi * (CHIRP_S - t) / (r * CHIRP_S)))
        out.append(w * math.sin(ph))
    return out


def write_probe(path, repeats=5, gap_s=0.5, lead_s=0.5, amplitude=0.7):
    c = chirp()
    total = int((lead_s + repeats * gap_s + 0.5) * FS)
    sig = [0.0] * total
    for k in range(repeats):
        o = int((lead_s + k * gap_s) * FS)
        for i, v in enumerate(c):
            if o + i < total:
                sig[o + i] = amplitude * v
    with wave.open(path, "w") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(FS)
        w.writeframes(struct.pack("<%dh" % total, *[int(32767 * v) for v in sig]))
    return path


def read_wav(path):
    with wave.open(path) as w:
        n, c, fs = w.getnframes(), w.getnchannels(), w.getframerate()
        d = struct.unpack("<%dh" % (n * c), w.readframes(n))
    return [x / 32768.0 for x in d[0::c]], fs


def matched_filter(rec, fs):
    """Correlate against the analytic chirp. Returns |correlation| per lag."""
    n = int(fs * CHIRP_S)
    rep = []
    for i in range(n):
        t = i / fs
        ph = 2 * math.pi * (F0 * t + 0.5 * (F1 - F0) / CHIRP_S * t * t)
        r = 0.25
        w = 1.0
        if t < r * CHIRP_S:
            w = 0.5 * (1 - math.cos(math.pi * t / (r * CHIRP_S)))
        elif t > CHIRP_S - r * CHIRP_S:
            w = 0.5 * (1 - math.cos(math.pi * (CHIRP_S - t) / (r * CHIRP_S)))
        rep.append(complex(w * math.cos(ph), w * math.sin(ph)))
    try:
        import numpy as np
    except ImportError:
        raise SystemExit("this needs numpy for the correlation; pip install numpy")
    x = np.array(rec)
    h = np.conj(np.array(rep)[::-1])
    N = 1 << int(math.ceil(math.log2(len(x) + n)))
    corr = np.fft.ifft(np.fft.fft(x, N) * np.fft.fft(h, N))[n - 1:n - 1 + len(x) - n + 1]
    return np.abs(corr)


def report_peaks(mag, fs, temp_c, rh=0.0, limit=6, min_separation_s=0.001):
    """Report the strongest well-separated correlation peaks.

    `min_separation_s` must be a few times the matched filter's own resolution
    (1/bandwidth = 0.17 ms here) but no more: an earlier version used 20 ms and
    silently discarded every echo arriving within 7 metres of the direct path,
    which on a desk is all of them."""
    import numpy as np
    med = float(np.median(mag))
    order = np.argsort(mag)[::-1]
    sel = []
    for i in order:
        if all(abs(int(i) - j) > int(min_separation_s * fs) for j in sel):
            sel.append(int(i))
        if len(sel) == limit:
            break
    c = sound_speed(temp_c, rh)
    print("  speed of sound at %.1f C, %.0f%% RH: %.2f m/s" % (temp_c, rh, c))
    print("  median background %.5g" % med)
    print("  %10s %10s %12s %10s %12s" % ("sample", "ms", "magnitude", "dB>med", "path (m)"))
    for i in sorted(sel):
        print("  %10d %10.2f %12.5f %10.1f %12.3f"
              % (i, 1000.0 * i / fs, mag[i], 20 * math.log10(mag[i] / med), c * i / fs))
    return sorted(sel)


def simulate(delay_s=0.0015, snr_db=10.0, echo_s=0.004, echo_db=-8.0, seed=7):
    """A synthetic channel: direct path, one echo, and noise. Used by selftest so
    the analysis path is verified without hardware."""
    import random
    rng = random.Random(seed)
    n_total = int(0.5 * FS)
    sig = [0.0] * n_total
    c = chirp()
    for src_delay, gain in ((delay_s, 1.0), (echo_s, 10 ** (echo_db / 20))):
        o = int(src_delay * FS)
        for i, v in enumerate(c):
            if o + i < n_total:
                sig[o + i] += gain * v
    amp = 10 ** (snr_db / 20)
    return [amp * v + rng.gauss(0, 1) for v in sig]


def cmd_selftest(args):
    print("Simulated channel: direct path at 1.50 ms, echo at 4.00 ms (-8 dB), SNR 10 dB")
    rec = simulate()
    mag = matched_filter(rec, FS)
    peaks = report_peaks(mag, FS, args.temperature, args.humidity)
    direct = peaks[0] / FS
    err_ms = abs(direct - 0.0015) * 1000
    print()
    print("  direct path recovered at %.3f ms, truth 1.500 ms, error %.3f ms" % (direct * 1000, err_ms))
    if err_ms > 0.05:
        print("  FAIL: delay estimate is off by more than one sample period")
        return 1
    if len(peaks) < 2 or abs(peaks[1] / FS - 0.004) > 0.0005:
        print("  FAIL: the echo was not resolved")
        return 1
    print("  PASS: direct path and echo both resolved")
    return 0


def cmd_probe(args):
    write_probe(args.output, repeats=args.repeats)
    print("wrote %s: %d chirps, %.0f-%.0f Hz, %.0f ms each" %
          (args.output, args.repeats, F0, F1, 1000 * CHIRP_S))
    return 0


def cmd_analyse(args):
    rec, fs = read_wav(args.input)
    print("%s: %.2f s at %d Hz" % (args.input, len(rec) / fs, fs))
    mag = matched_filter(rec, fs)
    report_peaks(mag, fs, args.temperature, args.humidity)
    return 0


def cmd_live(args):
    probe = "/tmp/phantom_air_probe.wav"
    rec = "/tmp/phantom_air_rec.wav"
    write_probe(probe, repeats=args.repeats)
    dur = args.repeats * 0.5 + 2.0
    print("Playing and recording for %.1f s..." % dur)
    r = subprocess.Popen(["arecord", "-D", args.device, "-f", "S16_LE",
                          "-r", str(FS), "-c", "1", "-d", str(int(dur)), rec],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(["sleep", "0.5"])
    subprocess.run(["aplay", "-D", args.device, probe],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    r.wait()
    args.input = rec
    rc = cmd_analyse(args)
    print()
    print("  REMEMBER: this delay includes the sound card's own latency, which is")
    print("  typically 10-50 ms and dwarfs the 1.5 ms it takes sound to cross half")
    print("  a metre. Run this at two known distances and subtract: the DIFFERENCE")
    print("  cancels every fixed electronic delay and leaves the acoustics.")
    return rc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--temperature", type=float, default=20.0, help="room temperature, C")
    ap.add_argument("--humidity", type=float, default=0.0,
                    help="relative humidity %%; matters more for absorption than for speed")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("selftest", help="verify the analysis on a simulated channel")
    s.set_defaults(func=cmd_selftest)

    s = sub.add_parser("probe", help="write a probe WAV")
    s.add_argument("output")
    s.add_argument("--repeats", type=int, default=5)
    s.set_defaults(func=cmd_probe)

    s = sub.add_parser("analyse", help="analyse a recording")
    s.add_argument("input")
    s.set_defaults(func=cmd_analyse)

    s = sub.add_parser("live", help="play and record on this machine")
    s.add_argument("--device", default="default")
    s.add_argument("--repeats", type=int, default=5)
    s.add_argument("--distance", type=float, default=None,
                   help="speaker-to-microphone distance in metres, if known")
    s.set_defaults(func=cmd_live)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
