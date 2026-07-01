## 1. The ontology, what things are

Everything reduces to four object types.

A **feature** is a name drawn from a universal, finite, innate set `F = {F₁, F₂, ..., Fₙ}`. Features are substance-free, they are labels with no phonetic content inside the phonology, and phonetic interpretation happens later in a separate transduction. `F` is the same for every language.

A **value** is drawn from `{+, −}`. There are exactly two. There is no third value.

A **valued feature** is a pair from `{+, −} × F`, written `+Voice` or `−Back`. Two valued features are **opposites** when they share a feature and differ in value, so `+F` and `−F` are opposites.

A **segment** is a set of valued features subject to one condition, consistency. A set `σ` is consistent iff there is no feature `F` such that both `+F ∈ σ` and `−F ∈ σ`. A segment simply cannot contain a feature and its opposite. This is the only well-formedness condition on segments.

Two formulations of a segment are isomorphic, and you should hold both in mind. The **set formulation** treats `/u/` as `{+Syllabic, +High, +Back}`. The **function formulation** treats a segment as a total function `σ : F → {∅, +, −}`, where `σ(F) = ∅` means underspecified for `F`. The two are interconvertible, the set is just the pairs where the function returns a value. Use the set formulation for stating rules and the function formulation when you reason about the full segment space.

The **segment space** is `Ω(F) = {σ | σ : F → {∅, +, −}}`, and its size is `|Ω(F)| = 3^|F|`, since each feature independently takes one of three states, plus, minus, or absent. That factor of three, not two, is the entire source of LP's analytic power. Underspecification is not an add-on, it is what you get for free the moment you allow the empty value.

## 2. The critical asymmetry, why absence cannot be targeted

This is the single most important conceptual point, and everything downstream depends on it.

A **fully specified** segment carries a value for every feature. An **underspecified** segment lacks a value for at least one. By convention underspecified segments are written with capitals, `/D/` for a coronal stop with no Voice value, `/T/` for one with no Strident value, and so on.

The rule is this. **You cannot write a rule that targets the absence of a value.** Absence is not a property, it is the lack of one, and there is no `∅Voice` you can name in a structural description. This is not a stipulation, it falls out of the natural-class definition in section 4, but state it to yourself now, because the whole framework is organized around obeying it. When a rule appears to single out an underspecified segment, it never does so directly. It targets a broad class that includes the underspecified segment and its fully specified relatives, and the operation itself sorts them out. The apparent targeting of `/D/` is always an illusion produced by the operator, never a reference to `/D/`.

## 3. Natural classes, the objects rules refer to

A **natural class** is a set of segments defined intensionally by a consistent feature set. Given a feature specification `C ⊆ ({+,−} × F)`, the class is

```
N(C) = { x ∈ Ω(F) | x ⊇ C }
```

A segment is in the class iff its feature set is a superset of `C`. This is the **subsumption** relation, and it is the whole story of class membership. Note the direction carefully, the class is defined by the *smaller* feature set, and it contains the *larger*, more specified segments.

**Specificity** is the immediate consequence. More features means a smaller class.

```
C₁ ⊇ C₂  implies  N(C₁) ⊆ N(C₂)
```

Adding a feature to a specification shrinks the class it picks out. The empty specification gives the universal class, `N(∅) = Ω(F)`, every segment. This matters later, since "adjacency" in SEARCH is exactly `N(∅)` as a terminator.

Three consequences you must internalize.

**Only a fully specified segment is a singleton natural class.** If `/t/` is fully specified, `N({t}) = {t}`, a class of one. But `/D/`, being underspecified, can never be a singleton class, because any specification that fits `/D/` also fits the more specified `/t/` and `/d/` that subsume it. This is the formal reason absence cannot be targeted. There is no consistent feature set whose superset-closure is exactly `{D}`.

**Most sets of segments are not natural classes.** With a two-feature system there are 9 segments, 2⁹ = 512 subsets, but only 9 natural classes. The other 503 cannot appear in any rule. The theory is precise about what a possible rule can mention, and this precision is the point.

**The Phonological Sawzall**, Reiss's name for the specificity corollary. If `f ⊂ g` as feature sets, there is no natural class containing the segment `f` while excluding the segment `g`. You cannot carve `g` out from under `f`. This is why you cannot target the less specified member to the exclusion of the more specified one, and it forces the underspecification analyses that are LP's signature.

## 4. The bracket convention, non-negotiable

LP takes brackets as type markers, and mixing them is a type error.

**Square brackets** `[ ]` denote a natural class, a set of segments, defined intensionally.

```
[+Syllabic, +High]  =  N({+Syllabic, +High})  =  the class of segments subsuming that set
```

**Curly braces** `{ }` denote a set of valued features, a plain set, used for the structural change.

```
{+High}  =  the feature set to be added or removed
```

The traditional SPE rule wrote both the target and the change in square brackets, and LP treats that as conflating two different set-theoretic types. Targets and environments are sets of segments, `[ ]`. Changes are sets of features, `{ }`. Get this right in every rule box, it is the first thing a practitioner checks, and your paper already does it correctly.

## 5. The circle notation

Now the circle notation you asked about, since you have the pieces to understand exactly what it does.

The problem it solves is this. Sometimes you want to target "every segment that carries at least the full specification of some particular segment `σ`," where `σ` is a richly specified, perhaps derived, segment. Writing out `σ`'s whole feature bundle in the environment is verbose and obscures that you mean "this segment's class." The circle abbreviates it.

**Circle notation.** For a segment `σ`, the circled symbol `Ⓢ` inside square brackets denotes the natural class whose defining specification is `σ`'s entire feature set.

```
[ Ⓢ ]  =  N(σ)  =  { x ∈ Ω(F) | x ⊇ σ }
```

It is ordinary natural-class notation with the defining set supplied by naming a segment rather than by listing features. For a fully specified `σ`, `[ Ⓢ ]` is the singleton `{σ}`. For a partially specified one, it is the class of everything at least as specified.

The circle removes exactly one layer of curly braces. Where `[ {+Cor, −Ant, +Dist, +Strid, +Cont} ]` would be the class defined by that bundle, if that bundle is the full specification of your derived postalveolar `Ⓢ`, then `[ Ⓢ ]` says the same thing compactly. This is precisely how your bleed rule uses it, `[ Ⓢ ]` in the environment picks out the derived postalveolar by its whole bundle without you re-listing `−Anterior, +Distributed, +Strident, +Continuant, +Coronal`.

When to reach for it. Use the circle when the class you want is "everything subsuming this specific, usually derived, segment," and spelling the bundle out would be long or would hide that you mean that segment's class. Do not use it when a short feature list already names the class you want, since then the circle is just indirection. In your basic S-palatalization the target `[+Coronal, +Strident, +Continuant]` is already a clean short class, so no circle. In your bleed environment the derived postalveolar's full bundle is what you mean, so the circle earns its place.

## 6. The operators, the content of Ω

Rules are built from a small set of operations on feature sets, the inventory called `Ω`. The program's central empirical question is what `Ω` contains. There are two live proposals. I will give the mainstream two-operator inventory first, since that is what your paper uses, then the one-operator alternative.

### 6a. Subtraction

Subtraction removes valued features.

```
A \ B  =  { cF | cF ∈ A ∧ cF ∉ B }
```

Everything in `A` survives except the exact valued features also in `B`. It is set difference. It applies vacuously when `B` shares nothing with `A`.

A subtlety your bleed rule uses. To remove a feature "at whatever value it bears," you subtract both polarities. `A \ {+F, −F}` strips `F` from `A` regardless of which value `A` had. Equivalently this is **projection** onto the complement of `{F}`, restricting the segment to the features other than `F`. So when your revised bleed says "delink Anterior, Distributed, Strident at whatever value," it subtracts both polarities of each, or projects away those three features. Both descriptions name the same operation.

Subtraction is total, always defined, and non-commutative. It is the operation that produces **derived underspecification**, taking a fully specified segment and making it underspecified mid-derivation, which is exactly what feeds your default fill-in.

### 6b. Unification

Unification adds valued features but only when the result stays consistent.

```
A ⊔ B  =  A ∪ { cF | cF ∈ B ∧ −cF ∉ A }
```

In words, `A ⊔ B` is everything in `A`, plus every feature of `B` whose opposite is not already in `A`. If `B` carries a value that contradicts one `A` already has, that feature of `B` is dropped, not added. Under the strict definition unification is a **partial operation**, formally undefined when `A ∪ B` would be inconsistent. But the rules built from it are **total**, because when unification is undefined the rule maps the segment to itself, vacuous application. This totality is what lets rules compose into a phonology.

The three outcomes, which are the entire mechanism of inalterability, using a target open for `F` versus specified for it:

```
{ }      ⊔ {+F}  =  {+F}     feature filling, the open segment takes the value
{+F}     ⊔ {+F}  =  {+F}     vacuous, already has it, no change
{−F}     ⊔ {+F}  =  {−F}     failure, opposite present, value refused
```

Row one is the underspecified segment mutating. Row three is the prespecified segment being **inalterable by unification failure**. The rule targets a class containing all three kinds of segment, and the operator alone decides which change. You never referred to the underspecified segment. It only *seems* like the rule singled it out. This sentence is the beating heart of every LP analysis, including yours, and you should be able to say it in your sleep.

**Feature-filling** is unification alone. **Feature-changing** is subtraction then unification, `d → D → t`, delete the offending value, then fill the new one. The two-step decomposition is not decoration, it makes two empirical predictions LP is proud of. **No bypassing**, a feature-changing rule cannot skip an underspecified segment, because to change `/d/` to `/t/` you must pass through `/D/`, and any underspecified `/D/` already in the input rides along. **No polarity rules**, no single rule maps `+F` to `−F` and `−F` to `+F` simultaneously, because deletion neutralizes the two before insertion can tell them apart. The absence of these rule types cross-linguistically is offered as evidence for the decomposition.

### 6c. The one-operator alternative, priority union

Reiss 2021 shows you can replace both operators with one, priority union, adapted from unification-based grammar. Priority union takes a **strict** argument and a **defeasible** argument and unions them, but on conflict the strict argument wins.

```
The priority union of A and B is the union of A with those valued
features of B whose opposites are not in A.
```

Written with a directional underarrow pointing strict-to-defeasible. Feature-filling makes the target strict, `a ⊔⃗ b`, so the target's existing values win and only its open features get filled. Feature-changing makes the change strict, `a ⃖⊔ b`, so the change overrides the target. One operator, total, non-commutative, and the filling-versus-changing distinction becomes a matter of which argument is defeasible.

The two proposals are extensionally equivalent over filling and changing, but they diverge on predictions. Subtraction can produce derived surface underspecification, priority union cannot. Priority union expresses exchange rules trivially, `[αF] ⃖⊔ {−αF}`, subtraction cannot. So the choice is empirical, not cosmetic. Your paper uses the subtraction-and-unification inventory with `⊔` and `\`, which is the mainstream choice, and your bleed rule's reliance on genuine subtraction is a live commitment to that inventory over priority union.

## 7. The third operator, segment insertion and deletion

Feature operations change segments' insides. A separate operator handles whole segments, using a null symbol `ε` and a mapping arrow `↦`.

**Insertion** links `ε` to a feature bundle and introduces a segment. **Deletion** maps a segment matching a natural class to `ε`.

```
Insertion:   ε ↦ {...} / ...
Deletion:    [ ... ] ↦ ε / ...
```

A crucial asymmetry, and a nice illustration of the whole system. Deletion targets are natural classes, so they obey specificity. **Delete the Rich**, if some but not all `X` delete in a context, the ones that delete must be *more* specified than the ones that stay, because only a more specified class can be carved out. This runs against the common intuition that weaker, less specified segments delete more easily, and LP embraces the reversal as a prediction. Note also that deleting a segment is not the same as removing all its features, since removing all features leaves the empty set, which is still a set, not nothing. `ε` is nothing. `{ }` is a segment with no specifications.

## 8. Rules, and how a grammar is a function

A **rule** has the schema

```
target  ω  change  /  environment
```

The target is a natural class in `[ ]`. `ω` is an operator from `Ω`, either `⊔`, `\`, or `↦`. The change is a feature set in `{ }`, or a bundle for insertion. The environment is stated with natural classes and a gap `__`. Reiss's slogan, **whatness is independent of whereness**, the features changed and the features of the environment need not overlap, so a rule may be "computationally arbitrary" and still well-formed, since phonetic naturalness lives in diachrony, not synchronic grammar.

Each rule is a **total function** from strings to strings, mapping each input segment to an output segment, identity where the rule does not apply. Because each rule is total, rules **compose**. The phonology of a language is the ordered composition of its rules,

```
phonology = rₙ ∘ ... ∘ r₂ ∘ r₁
```

Ordering is extrinsic and does real work. Your bleed-before-assibilation ordering is exactly this, bleed must delink Strident before assibilation looks for a target, or the wrong segment mutates.

## 9. SEARCH, the environment as an algorithm

The operators above are the *what*. SEARCH is the *where*. It is a separate module, and the two are explicitly complementary, the operator papers abstract away from environments and the SEARCH papers abstract away from intrasegmental change. Azadegan 2026 is the paper that welds them.

A SEARCH rule is a quintuple `⟨Inr, Trm, Dir, Out, Cnd⟩`.

**Inr**, the initiator class, the targets, the segments that launch a scan. **Trm**, the terminator class, the segments that can halt a scan. **Dir**, the direction, `Left` or `Right`. **Out**, the output function, what change applies, built from the operators of section 6. **Cnd**, a condition, a natural class the terminator must satisfy for the change to fire.

The mechanism. Each initiator launches a scan in direction `Dir`. The scan runs to the nearest segment in `Trm`. That terminator is then checked against `Cnd`. If it satisfies `Cnd`, `Out` applies to the initiator. If not, nothing happens, and crucially the scan does *not* continue past the terminator, it has already stopped.

Formally, for position `i` searching right, `Θ_R(i)` is the least `j > i` with `s_j ∈ Trm`, undefined if none exists. Position `i` is licensed iff `s_i ∈ Inr`, `Θ_R(i)` is defined, and `s_{Θ_R(i)} ∈ Cnd`. Licensed positions get `Out`, all others are copied unchanged, and all outputs are computed simultaneously from the input, so within one rule nothing feeds or bleeds.

The two design payoffs, both of which your paper uses.

**Adjacency is `Trm = N(∅)`.** Make the terminator the universal class and the scan halts at the immediate neighbor every time. Locality is not a separate device, it is the degenerate case where the terminator is maximally general. This is the "adjacency is opaqueness" thesis.

**Transparency and blocking come from where the discriminating class sits.** Put the discriminating features in `Trm`, and non-members are simply not terminators, so the scan passes through them, they are *transparent*. Put the discriminating features in `Cnd` with a broad `Trm`, and every segment halts the scan, so a segment that fails `Cnd` *blocks*. Same machinery, opposite effect, selected by the terminator-versus-condition choice. Your S-palatalization puts the postalveolar class in `Trm`, so the intervening stop is transparent and `/s/` reaches the vowel across it. Your dorsal palatalization puts the front-vowel condition in `Cnd` with a next-segment `Trm`, so an intervening stop halts the scan and blocks. One parameter, terminator breadth, generates both.

The initiator-first design has a formal virtue worth knowing. Mapping each initiator to its single terminator is a function, since each initiator has exactly one nearest terminator. Mapping a trigger to its targets is not a function, since one trigger can serve many targets. So "recipient as initiator" buys a well-defined mapping, which is why the framework insists the target does the looking, not the trigger.

## 10. The complete toolkit, assembled

To implement an LP analysis you need, in order.

A universal feature set `F`, and lexical entries as consistent, possibly underspecified feature sets, capitals for the underspecified ones. A decision about which segments are prespecified and which underspecified, driven by the three-way-behavior diagnostic, if two surface segments show three behaviors, the third is underspecified. Natural classes for every target and environment, written `[ ]`, using the circle `[ Ⓢ ]` where you mean "everything subsuming this segment." Changes written `{ }`. Operators from `Ω`, unification for filling, subtraction for delinking, the two composed for feature-changing, segment `↦ ε` for deletion. Default fill-in rules, late and context-free, to supply unmarked values to whatever the mutation rules left open, the step where your `/D/` spirantizes and your unmutated stops surface plain. An extrinsic rule ordering. And, for the environments, either simple `__` gaps read as SEARCH terminators, or explicit quintuples when the locality is doing real work.

The recurring analytic move, which is the reason to hold all of this, is one shape. A pattern looks like it needs morphology, diacritics, or exceptions. You find three behaviors in two segments, posit underspecification, write one broad-target unification rule, and the prespecified segments fall out as inalterable by unification failure while the underspecified ones mutate. The exceptionality was never in the grammar. It was in the lexicon. Your Romanian `/k/` versus `/K/`, `/t/` versus `/T/`, `/s/` versus `/S/` are each one instance of exactly that move, and P collects the three inalterability conditions into a single statement that also drives the allomorphy. That is LP, and your paper is a clean, multi-place-of-articulation instance of it.

Here is the same system explained as conceptual machinery, what each piece is *doing* and why the pieces fit together the way they do.

## The founding decision, three states instead of two

Every phonological theory decomposes segments into features. LP's one consequential move is to let a feature be in three states rather than two, present-and-positive, present-and-negative, or absent altogether. That third state, absence, is the entire engine. It is not a diacritic bolted on to handle exceptions, it is what you automatically get once you stop requiring segments to answer every feature question. A theory that fills in every feature is making an extra stipulation, completeness, and LP simply declines to make it.

The reason this matters is that absence gives you a way to encode a difference between two segments that behave differently without inventing a feature to distinguish them. Two segments can be identical in every specified feature and differ only in that one of them has a value where the other has a hole. That hole is invisible to the ear, it surfaces as some default, but it is visible to the grammar, because rules interact with holes differently than with filled values. So underspecification is a channel for lexical information that does not show up directly on the surface but changes how a segment participates in processes. That channel is what LP spends all its power exploiting.

## Why you can never point at a hole

The deepest structural fact in the system is that a rule cannot mention the absence of a value. You can target segments that have `+F`, you can target segments that have `−F`, but you cannot target segments that lack `F`, because "lacking `F`" is not a property, it is the failure to have one, and the machinery for describing rule targets only knows how to describe things by the properties they possess.

This is not a rule someone imposed. It falls out of how targets are described. A rule describes its target by listing some features, and the target is everything that has *at least* those features. Listing features can only ever enlarge the description downward, toward more-specified segments, never carve out a less-specified one. If you try to describe "the segment missing a Voice value," any feature list you write that fits it also fits the fuller segments that have Voice filled in, because those fuller segments have everything your list mentions plus more. So the underspecified segment can never be isolated by description. It always drags its more-specified relatives into the target with it.

The immediate corollary is that a bare underspecified segment is never a class of its own. A fully specified segment can be, you can describe it so tightly that only it qualifies, but an underspecified one cannot be pinned down, because there is always something more specified sitting on top of it that your description also catches. This one asymmetry, fully specified segments can be singled out and underspecified ones cannot, is the pressure that shapes every LP analysis. It forces you to write rules that aim broadly and then let something else do the sorting.

## The sorting is done by the operation, not the target

Since a rule cannot aim narrowly at the underspecified segment, LP relocates the selectivity. The rule aims at a broad class that contains both the segment you want to affect and the segments you want to spare, and then the *operation* the rule performs succeeds on some members and fails on others. The rule does not decide who changes. The arithmetic of the change decides.

The operation that does this is unification, and its whole job is to add information only where there is room for it. Unification tries to write a value into a segment. If the slot is empty, the value goes in. If the slot already holds the same value, nothing happens, harmlessly. If the slot already holds the opposite value, the write is refused, because you would create a contradiction, a segment that is both plus and minus the same feature, which is not a legal segment. That refusal is not an error to be handled, it is the selectivity you wanted. The segments that get changed are exactly the ones with an empty slot. The segments that resist are exactly the ones already committed to the opposite value. The rule reached for all of them equally, and the operation quietly spared the committed ones.

This is why LP keeps saying it only *seems* like the rule targets the underspecified segment. The targeting is an illusion produced downstream, by unification failing on the specified segments. Inalterability, a segment's immunity to a process, is not a property the segment is marked with, and it is not a condition the rule checks. It is just what happens when you try to write a value into a slot that already contradicts it. A prespecified segment is inalterable because there is no room to change it, and that is the entire explanation.

Once you see this, the analytic recipe writes itself. If two surface segments show three distinct behaviors, one of them must secretly be two things, a committed version and an uncommitted version. The committed version carries the value and resists the rule. The uncommitted version lacks it and undergoes. You do not need morphological features, exception diacritics, or special rule marking. You need one lexical distinction, filled versus empty, and one broad rule whose operation sorts them. The apparent exceptions were never grammatical, they were lexical all along.

## Building change out of two moves

Unification only adds, it never overwrites, so on its own it can only fill holes. Real languages also flip values, a voiced thing becomes voiceless. LP refuses to give itself a single operation that flips a value, and instead builds flipping out of two moves, first erase, then fill. Erasing is a separate operation, subtraction, which removes a value and leaves a hole. Filling is unification again. So changing `+F` to `−F` means first subtracting the `+F`, which makes the segment temporarily underspecified, then unifying in `−F`, which the now-empty slot accepts.

This looks like a roundabout way to do something simple, but the roundaboutness makes predictions, and the predictions are the payoff. Because a value-flip must pass through an underspecified stage, the flipped segment becomes momentarily identical to a segment that was underspecified to begin with. At that moment the grammar cannot tell them apart, and anything that happens next happens to both. This means a flipping rule can never quietly skip over an underspecified segment to get at a specified one, because the moment it operates it has erased the very distinction that would let it skip. And it means you cannot write a rule that flips plus-to-minus and minus-to-plus at the same time, because the erasure step wipes out which one each segment started as before the filling step could tell them apart and give them opposite new values. Languages do seem to lack exactly these rule types, and LP offers their absence as evidence that change really is erase-then-fill rather than a primitive flip.

There is a competing account that does it all with one operation, priority union, where you designate one side as the winner in a conflict. If the target wins, you get filling, existing values are protected and only holes are filled. If the change wins, you get flipping, the change overrides whatever was there. This is tidier, one operation instead of two, but it is not free, it makes different predictions about whether a segment can be left underspecified after a process and about whether value-swapping rules can exist. So the choice between the two-operation and one-operation stories is a real empirical bet, not a matter of taste. Committing to genuine erasure, as an analysis does whenever it needs a segment to become underspecified in the middle of a derivation, is a commitment to the two-operation inventory.

## Defaults clean up afterward

Because rules leave holes, and because erasure creates new holes mid-derivation, something has to fill the holes before the segment reaches the surface, since the surface has no room for underspecified segments, they have to be pronounced as something. That job goes to late, unconditional fill-in rules that supply the unmarked value to any slot still empty. These run after the interesting rules and simply mop up.

This division of labor is doing conceptual work, not just tidiness. It means a segment's surface form can be assembled from two sources, whatever the active rules put in, plus whatever the defaults supply for what is left. Two segments that end a derivation with different holes will get different default fillings and surface differently, even though no rule ever explicitly distinguished them. This is how one underlying contrast, a hole here versus a value there, can cash out as a surface contrast that looks like it needed its own rule but did not. The contrast rides through the derivation as a difference in what is left open, and the defaults turn that difference into audible output at the end.

## Separating what changes from where it changes

So far this is all about the change itself. The other half of a rule is the environment, the condition on where the change happens, and LP insists these two halves are logically independent. What a rule does to a segment need not have anything featurally in common with what triggers it. A rule can raise a vowel before a nasal without copying anything from the nasal, and the theory treats "raise before a nasal" and "nasalize before a nasal" as the same kind of object, a change conditioned by an environment, with no special status for the case where the change happens to echo the trigger. Assimilation, where the change does echo the environment, is not a privileged category, it is just the accidental case where the what and the where share features. The framework refuses to build assimilation in as a basic notion, because doing so would treat a coincidence as a mechanism.

## The environment as a search, and locality as a side effect

The environment is where LP makes its most distinctive move. Instead of stating environments as static templates, it states them as a search. A target segment launches a scan in some direction, the scan travels until it hits a segment of a designated stopping-class, and then the grammar inspects what it hit. The change fires or does not depending on what the scan landed on.

The reason to do it this way is that it dissolves locality as a separate concept. In most theories, adjacency and action-at-a-distance are two different mechanisms, and you need extra machinery to say which segments an intervening thing is allowed to block. In the search picture, there is only one mechanism, and locality is a setting inside it. If you declare that *any* segment can stop the scan, the scan always stops at the immediate neighbor, and you have adjacency, without adjacency being its own device. If you declare that only certain segments stop the scan, the scan sails past everything else until it finds one, and you have action-at-a-distance. Adjacency is not the basic case with distance as an exotic elaboration, it is the other way around, distance is the default and adjacency is what you get by making the stopping-class maximally permissive.

This also dissolves transparency and blocking, the two ways an intervening segment can behave. Whether something in the middle is invisible to a process or halts it is not a property of that middle segment, it is a consequence of how you set up the search. There are two slots where selectivity can live. You can put it in the stopping-class, so that only the right kind of segment ever halts the scan, in which case everything else is passed over silently and counts as transparent. Or you can make the stopping-class permissive so the scan halts at the first thing it meets, and put the selectivity in a separate check applied *after* stopping, in which case a segment that halts the scan but fails the check ends the search unsuccessfully and counts as a blocker. Same intervening segment, opposite behavior, and the difference is entirely in whether the discriminating property was placed in the stopping-condition or in the post-stop check. Transparency and opacity are not things segments have, they are shadows cast by where you located the selectivity.

There is a reason the target does the searching rather than the trigger reaching out to its targets, and it is not arbitrary. If each target looks for its own single stopping-point, the relationship is clean, every target has exactly one thing it found. If instead a trigger reached out to all the targets it could affect, one trigger could be responsible for many targets, and the relationship stops being a tidy one-to-one map. Making the target the active party keeps the geometry of the rule well-behaved, each target paired with the one segment its search terminated on.

## How the whole thing coheres

Step back and the pieces are all serving one idea. Underspecification gives you a hidden lexical channel, a difference between a hole and a filled value that the ear cannot detect but the grammar can. The impossibility of pointing at a hole forces rules to aim broadly. Unification, by adding only into holes and refusing to overwrite, turns that broad aim into precise selectivity, mutating the uncommitted segments and sparing the committed ones, with inalterability falling out for free as failed writes rather than as marked exceptions. Subtraction lets you manufacture holes on purpose, which both enables value-flipping as erase-then-fill and predicts the rule types that never occur. Defaults convert leftover holes into surface pronunciations, letting one underlying hole-versus-value contrast surface as an audible difference no single rule created. The independence of what from where keeps the change and its trigger from being confused for one mechanism. And the search-based environment makes locality, transparency, and blocking into settings of one procedure rather than separate devices, so that whether a segment intervenes harmlessly or fatally is just a matter of where you parked the selectivity.

The through-line is that LP keeps taking things other theories treat as primitives, exception-marking, feature-filling versus feature-changing as distinct rule types, adjacency, transparency, blocking, and showing that each is a downstream consequence of a smaller set of moves, three feature-states, add-only unification, erasure, and directional search. The payoff for your Romanian is that a wall of apparent morphological conditioning and lexical exceptions reduces to one lever, whether a given obstruent carries its distinguishing value or leaves it open, pulled once per place of articulation, with everything else, the inalterability, the allomorph selection, the locality of the cluster effects, following from the machinery rather than being stipulated.