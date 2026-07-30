#!/bin/bash
# End-to-end GUI verification of the Phased Light Curve widget.
set -u
cd /home/krutkowski/Pulpit/ihsnpeaks-dev/ihsnpeaks/downstream/lc-qt || exit 1

./build/lc-qt ../../test_data/OGLE-BLAP-035.dat >/tmp/e2e_app.log 2>&1 &
APP=$!
sleep 3
WID=$(xwininfo -root -tree | grep '"lc-qt v' | awk '{print $1}')
echo "WID=$WID PID=$APP"
if [ -z "$WID" ]; then echo "NO WINDOW"; cat /tmp/e2e_app.log; exit 1; fi

alive() { kill -0 "$APP" 2>/dev/null && echo "alive" || echo "DEAD"; }

# count saturated-blue pixels in a region of a screenshot (hover highlight / scatter)
blue_in() { # file x y w h
  convert "$1" -crop "$4x$5+$2+$3" txt:- 2>/dev/null | tail -n +2 | \
    awk '{split($2,a,","); r=a[1]+0; g=a[2]+0; b=a[3]+0; if (b>140 && b-r>70 && b-g>30) c++} END{print c+0}'
}

echo "== step 1: locate Period Search button row via hover highlight =="
FOUND_Y=""
for Y in $(seq 500 8 720); do
  xdotool mousemove --window "$WID" 115 "$Y" 2>/dev/null
  sleep 0.2
  import -window "$WID" /tmp/hover.png 2>/dev/null
  C=$(blue_in /tmp/hover.png 20 $((Y-16)) 400 32)
  if [ "$C" -gt 300 ]; then FOUND_Y=$Y; echo "button row found at Y=$Y (blue=$C)"; break; fi
done
if [ -z "$FOUND_Y" ]; then echo "BUTTON ROW NOT FOUND; app=$(alive)"; exit 1; fi

echo "== step 2: scan X at the row to find the 2nd button (IHS) =="
RUNS=""
STATE=0
for X in $(seq 14 4 430); do
  xdotool mousemove --window "$WID" "$X" "$FOUND_Y" 2>/dev/null
  sleep 0.12
  import -window "$WID" /tmp/hover.png 2>/dev/null
  C=$(blue_in /tmp/hover.png $((X-30)) $((FOUND_Y-16)) 60 32)
  if [ "$C" -gt 300 ]; then
    if [ "$STATE" -eq 0 ]; then RUN_START=$X; STATE=1; fi
    RUN_END=$X
  else
    if [ "$STATE" -eq 1 ]; then RUNS="$RUNS $(( (RUN_START+RUN_END)/2 )),"; STATE=0; fi
  fi
done
echo "button centers:$RUNS"
IHS_X=$(echo "$RUNS" | tr ',' ' ' | awk '{print $2}')
if [ -z "$IHS_X" ]; then echo "IHS BUTTON NOT FOUND"; exit 1; fi
echo "IHS center at ($IHS_X, $FOUND_Y)"

echo "== step 3: click IHS and wait for the periodogram =="
xdotool mousemove --window "$WID" "$IHS_X" "$FOUND_Y" 2>/dev/null
sleep 0.2
xdotool click --window "$WID" 1 2>/dev/null
sleep 20
import -window "$WID" /tmp/e2e_spectrum.png 2>/dev/null
S=$(blue_in /tmp/e2e_spectrum.png 60 640 380 190)
echo "spectrum region blue after IHS: $S (app=$(alive))"

echo "== step 4: click the spectrum plot to set the pivot =="
SPEC_Y=$(( (FOUND_Y + 40 + 838) / 2 ))
xdotool mousemove --window "$WID" 226 "$SPEC_Y" 2>/dev/null
sleep 0.2
xdotool click --window "$WID" 1 2>/dev/null
sleep 1
import -window "$WID" /tmp/e2e_phased.png 2>/dev/null
P=$(blue_in /tmp/e2e_phased.png 470 40 410 410)
echo "phased plot blue after pivot click: $P (app=$(alive))"

echo "== step 5: move the pivot (click another spectrum frequency) =="
xdotool mousemove --window "$WID" 120 "$SPEC_Y" 2>/dev/null
sleep 0.2
xdotool click --window "$WID" 1 2>/dev/null
sleep 1
import -window "$WID" /tmp/e2e_phased2.png 2>/dev/null
P2=$(blue_in /tmp/e2e_phased2.png 470 40 410 410)
DIFF=$(compare -metric AE /tmp/e2e_phased.png /tmp/e2e_phased2.png null: 2>&1)
echo "phased plot blue after 2nd pivot: $P2; differing pixels: $DIFF (app=$(alive))"

kill "$APP" 2>/dev/null
echo "== done =="
exit 0
