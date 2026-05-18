# Stimulus Alignment Comparison: Your Script vs Master Script

## CRITICAL DIFFERENCES FOUND

### 1. **STIMULUS END TIME CALCULATION** ⚠️ MAJOR ISSUE
**Your Proven Method:**
```matlab
stim_total_time = data.Stimuli(2).stimulus_trial_t * data.Stimuli(2).trials
stim_start = data.Stimuli(2).TimeStimulusFrame(1)
stim_end = stim_start + stim_total_time
```
- Uses `TimeStimulusFrame(1)` as start
- Calculates end: `start + (stimulus_trial_t * num_trials)`
- Explicit duration calculation

**Master Script Current:**
```matlab
stim_time_values = Stimuli(i).TimeStimulusFrame;
time_start = min(stim_time_values(:));
time_end = max(stim_time_values(:));
```
- Uses min/max of entire TimeStimulusFrame array
- **Could give different results!**
- May include overhead/buffer times

**ISSUE:** These can calculate different time windows. Your method is more reliable because it's explicit about duration.

---

### 2. **BOUNDARY CONDITIONS** ⚠️ MEDIUM ISSUE
**Your Method:**
```matlab
Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
```
- Uses strict `>` and `<` (excludes boundary frames)

**Master Script (Current):**
```matlab
frame_idx = find(time_vector >= time_start & time_vector <= time_end);
```
- Uses `>=` and `<=` (includes boundary frames)

**ISSUE:** Your method excludes boundary frames; Master includes them. This affects which calcium frames are captured.

---

### 3. **TRIAL CHUNKING LOGIC** ⚠️ MEDIUM ISSUE
**Your Method:**
```matlab
num_trials = data.Stimuli.trials
chunk_size = length(data.Stimuli(1).TimeStimulusFrame)/num_trials;

for idx = 1:no_stimuli 
    for no_trial = 1:num_trials
        desired_vector = data.Stimuli(idx).TimeStimulusFrame;
        desired_vector = desired_vector(chunk_size * (no_trial -1) + 1 : chunk_size * no_trial);
        data.Stimuli(idx).TrialTimes{no_trial} = desired_vector;
    end
end
```
- Divides TimeStimulusFrame array into chunks (one per trial)
- Creates TrialTimes cell array for trial-by-trial access
- Used for subsequent trial-level analysis

**Master Script Current:**
- Does NOT chunk trials
- Processes entire stimulus duration as one unit
- No trial-level organization

**ISSUE:** Master script loses trial structure. You later use this for trial-by-trial averaging and correlation analysis.

---

### 4. **TIMECA ROW USAGE**
**Your Method:**
- Consistently uses `TimeCa(1, :)` for time vector
- Sometimes uses `TimeCa(2, :)` for different purposes
- Both rows available

**Master Script:**
- Uses `TimeCa(1, :)` (correct)
- Handles both 1D and 2D formats

**STATUS:** ✅ Aligned

---

## RECOMMENDATION

The Master script needs these fixes:

1. **Calculate stim_end explicitly** instead of using min/max
2. **Change back to strict `>` and `<`** (exclude boundaries)
3. **Add trial chunking logic** to create TrialTimes
4. **Preserve trial structure** for downstream trial-level analysis

These changes would make Master_Stimulus_Analysis match your proven, working approach.
