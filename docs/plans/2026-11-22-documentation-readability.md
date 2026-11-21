# Current documentation readability review

Reviewed104 current Markdown documents:root README,runtime guide,SALMONpy human/AI guides,and example indexes. Historical plans and frozen source documents remain evidence. Rewrote onboarding,coupling,API,saved-state and installation introductions around user goals and concrete execution methods. Removed internal naming from explanatory prose while retaining public callable names,paths,CLI flags and schema identifiers. Removed human-facing topic/revision bookkeeping,corrected the outdated API-version introduction,and kept developer linking details after the user workflow.

The main choice is now explicit:use the installed LD/SBE components;run a SALMON executable in a separate process;or call shared libraries one step at a time from Python. Failure/retry and shared-HSE introductions explain their purpose before API names. No numerical implementation changes were made in this review.

Checks:104 documents,zero unexplained native labels outside identifiers/paths,zero broken local links,zero unintended executable-example changes. Only comments and the documented current one-argument source-preparation command change in examples. AI guide bodies and plain-language search are synchronized. Manual11 and public-API5 tests pass. Numerical source inventory matches the previously validated reader runtime. No solver rebuild or CI was required.

Evidence:calc2026/salmonpy-documentation-readability-20261122/summary.json. Uncommitted.

## Follow-up:current status and navigation

The Japanese landing page now explains units,LD response,state reuse and execution environment in Japanese. Its obsolete per-step-state-reconstruction claim was corrected against the persistent state code,and the directional LD class is named explicitly. Paired coupling guides provide purpose-based links and distinguish create_state/open_model/open_registry in ordinary language. Superseded claims that Python shared models were unconnected were replaced with the current dedicated guides. Internal communication/resource and validation records are retained in developer details. First-use guidance and AI task searches use the same everyday descriptions;formal API IDs remain unchanged.

Fresh checks:manual11 and public API5 pass;all current guide links and executable-example comparisons pass. No numerical source change or additional CI/build was needed. Evidence:calc2026/salmonpy-documentation-readability-20261122/followup-summary.json. Uncommitted.
