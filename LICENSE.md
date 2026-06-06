# TAGMOS Public Research License v1.0

**Effective date**: 2026-06-06
**Licensor**: Wellmicro S.r.l. — Via Antonio Canova 30, 40128 Bologna, Italy
**Software**: TAGMOS framework v4.5.13 FORMULA A FULL16 — public research bundle

---

## 1. Grant of license

Subject to the terms and conditions of this Licence, Wellmicro S.r.l. ("Licensor")
hereby grants you ("Licensee") a worldwide, royalty-free, non-exclusive,
non-transferable, non-sublicensable licence to use, reproduce, modify and
prepare derivative works of the TAGMOS framework, schema, command-line
classifier, calibration recipe and example code (the "Software") contained
in the present public release bundle, **for the sole purpose of academic,
not-for-profit research**, subject to the conditions below.

## 2. Permitted use · academic non-commercial research only

You may:
- run the Software on metagenomic data of cohorts under your control;
- recalibrate the framework on your own cohort following the calibration
  recipe in `03_CALIBRATION_RECIPE/`;
- modify the Software for the internal purposes of your research group;
- prepare derivative works for the purpose of academic publication, public
  benchmark comparison, methodological extension or curriculum / teaching;
- publish numerical results, figures, supplementary tables and source code
  derived from the Software in peer-reviewed journals, preprint servers,
  conference proceedings, theses and other non-commercial academic venues.

## 3. Prohibited use · commercial activities

You may NOT, without an explicit separate commercial licence agreement
from Wellmicro S.r.l.:

a. use the Software, in whole or in part, in any product, service or
   activity intended for, or resulting in, commercial gain — including but
   not limited to: commercial diagnostic services, fee-for-service clinical
   reports, software-as-a-service platforms, microbiome-testing consumer
   products, paid clinical decision support, paid biomarker discovery
   pipelines, or any operation that bills end-users, patients, healthcare
   providers, insurance payers or pharmaceutical sponsors for outputs that
   incorporate the Software;

b. integrate the Software, in whole or in part, into a proprietary closed
   product or platform offered for sale, subscription, licensing fee, royalty
   or equivalent commercial consideration;

c. use the Software to generate evidence in support of regulatory
   submissions for commercial diagnostic devices, in vitro diagnostic
   software, software-as-a-medical-device (SaMD) clearances or any other
   regulatory pathway intended to enable commercial sale of a product
   incorporating the Software outputs;

d. publish or disclose any modification of the Software, derivative work
   or recalibrated instance under a licence that is permissive with respect
   to commercial use (e.g. MIT, BSD-2, BSD-3, Apache 2.0 or equivalent),
   unless explicit prior written consent of Wellmicro S.r.l. has been
   obtained.

If your intended use is commercial in any of the senses above, you are
required to obtain a separate commercial licence by contacting
`andrea.castagnetti@wellmicro.com`.

## 4. Mandatory attribution and citation

Any publication, preprint, presentation, dataset, software repository,
poster, abstract, thesis or other work — academic or otherwise — that uses
the Software, any modification of the Software, any output of the Software,
any recalibrated instance of the Software or any methodology described in
the Software bundle, **must cite the TAGMOS framework as follows**.

### Primary citation (mandatory)

> Soverini M., Lotfollahdzadeh A., di Rito L., Viciani E., Padella A.,
> Santacroce B., Marcante A., Monaldi C., Velichevskaya A., Castagnetti A.
> *Functional multi-axis decomposition of the human gut microbiome: an
> operational definition of eubiosis and dysbiosis.* bioRxiv 2026.
> doi: [TBD — populate with the bioRxiv DOI on publication]

### Recommended additional citation if specific extensions are used

If the user has used the cross-pipeline robustness battery, the
Wood-Ljungdahl ancestral-preservation analysis, the Crohn-vs-UC
sub-phenotype discrimination or the cross-population transferability
analysis, the corresponding Supplementary Note (SN1 through SN5) should
also be cited.

### Specific attribution language to include in Methods

Authors are required to include in the Methods section of the publication a
sentence equivalent to:

> *"Substrate-functional decomposition was computed using the TAGMOS
> framework (Soverini et al., bioRxiv 2026), recalibrated on the present
> cohort following the public calibration recipe (TAGMOS_PUBLIC_BUNDLE_v10_5)
> under the TAGMOS Public Research Licence v1.0."*

### Specific attribution language for derivative work

If you publish a *modification* of the framework (new axis, new EC entry,
re-engineered formula) you must state:

> *"This work extends/modifies the TAGMOS framework (Soverini et al.,
> bioRxiv 2026). The base architecture, EC dictionary and calibration recipe
> derive from the public release of TAGMOS v4.5.13. Modifications are
> described in §[…] of the present paper."*

Failure to provide the citation and attribution as specified above
constitutes a breach of this Licence and terminates your rights under
Section 1.

## 5. What this licence does NOT cover

This Licence specifically does NOT grant you any rights to:

- the proprietary frozen calibration parameters (per-axis mean, standard
  deviation and quartile / tertile cuts) of the TAGMOS Italian RWE
  calibration anchor — these are held as proprietary calibration parameters
  of Wellmicro S.r.l. and are not released in this bundle (see
  `03_CALIBRATION_RECIPE/README.md` for the rationale and for the calibration
  recipe that allows you to compute your own equivalent parameters on your
  own cohort);
- the underlying Wellmicro Metagenomics Pipeline (WMP) upstream profiler
  source code — only the framework, schema and bundled scripts are
  distributed under this Licence;
- the Wellmicro 6,508-subject Italian real-world evidence cohort
  subject-level data — only summary statistics published in the
  accompanying preprint are released;
- the Wellmicro trademark, brand assets, logo, font and visual identity —
  no permission is granted by this Licence to use these in association
  with any derivative work without a separate written authorisation.

## 6. Patent and intellectual-property reservation

The TAGMOS framework is the subject of one or more patent applications,
filed and held by Wellmicro S.r.l., covering both the framework architecture
and certain clinical-decision-support implementations. This Licence does
NOT grant any rights, express or implied, under those patents for any
commercial use as defined in Section 3. Academic, non-commercial research
use as defined in Section 2 is hereby authorised by the Licensor and does
not require additional patent licensing during the term of this Licence.

## 7. No warranty

The Software is provided "as is", without warranty of any kind, express
or implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and non-infringement. In no event shall
the Licensor be liable for any claim, damages or other liability arising
from, out of or in connection with the Software or the use or other dealings
in the Software.

## 8. Termination

This Licence is granted for as long as you comply with its terms. Failure
to comply with any term (in particular Section 3 prohibited commercial use,
or Section 4 mandatory citation) terminates your rights under this Licence
automatically. Upon termination, you must cease use of the Software and
destroy all copies in your possession.

## 9. Governing law

This Licence is governed by the laws of Italy. Any dispute arising under or
in connection with this Licence shall be subject to the exclusive
jurisdiction of the courts of Bologna, Italy.

## 10. Contact

For commercial licensing requests, patent licensing requests, or
clarifications on the scope of permitted academic use:

**Andrea Castagnetti** · Chief Scientific Officer · Wellmicro S.r.l.
email: andrea.castagnetti@wellmicro.com

---

*By using the Software you accept the terms of this Licence in full.*

*TAGMOS Public Research Licence v1.0 · 2026-06-06*
