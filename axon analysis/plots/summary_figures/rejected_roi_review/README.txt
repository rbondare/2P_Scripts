ROI review queue: for every aggregated session with rejected suite2p ROIs, a ranked thumbnail sheet of the rejected ROIs that look most like axons / boutons, and a decision sheet to fill in.
usage: python review_queue.py [session ...]          (default: every session that has a cache/rejroi_<session>.pkl table; build the tables with rejected_rois.py first)
output (plots/summary_figures/rejected_roi_review/):  review_<session>.pdf   labels_review_<session>.csv   review_queue_summary.csv   README.txt

Queue = rejected ROIs (iscell 0) that are not hidden constituents of a merged ROI (suite2p inmerge), do not overlap an accepted ROI, and are activity-like: they pass the axon transient test used for the
accepted ROIs OR have skewness > 0.5 (transient-rich trace).  Ranked by an activity score = mean percentile rank (within the session's rejected ROIs) of skewness, slow variance (transient test), local contrast
and F / Fneu.  At most MAX_Q ROIs per session.  The thumbnails show only morphology and the dF/F trace (no behaviour), so the decision cannot be driven by the coupling that is being tested.
Rejected ROIs with a median neuropil-corrected fluorescence <= 5 counts (dim) are not part of the tables and therefore not in the queue.
Fill in the 'decision' column of labels_review_<session>.csv with accept / reject / unsure (blank = no decision), then run apply_curation.py.

After filling in the decisions: python apply_curation.py [--apply] [session ...]   (dry run without --apply)
