#!/usr/bin/env python3
"""Templated arithmetic word problems with known integer answers.

Templates are structurally varied (multi-step, remainder, comparison, unit,
rate, division-with-leftover) so that persona reading styles produce
DIFFERENT failure modes rather than different overall accuracy — the
compatibility manipulation depends on that.

Each item: {"id", "template", "text", "answer"}. Deterministic under --seed.
"""
import argparse
import json
import random


def gen_items(seed: int, per_template: int) -> list[dict]:
    rng = random.Random(seed)
    items = []

    def add(template, text, answer):
        items.append({"id": f"{template}-{len(items):04d}", "template": template,
                      "text": text, "answer": int(answer)})

    for _ in range(per_template):
        # multi-step buy/lose (2-digit multiply)
        a, b, c = rng.randint(11, 48), rng.randint(7, 29), rng.randint(13, 59)
        add("boxes", f"A shop sells pencils in boxes of {a}. A school buys {b} boxes, then {c} pencils break and are thrown away. How many pencils remain? Answer with only the final integer.", a * b - c)
        # division with leftover, larger numbers
        n, k = rng.randint(151, 480), rng.randint(7, 19)
        add("seating", f"A hall puts {k} chairs per row. There are {n} guests. How many FULL rows can be seated? Answer with only the final integer.", n // k)
        # three-spend wallet
        start, s1, s2, s3 = rng.randint(180, 420), rng.randint(17, 68), rng.randint(23, 74), rng.randint(11, 49)
        add("wallet", f"Mara has ${start}. She spends ${s1} on lunch, ${s2} on a ticket, and ${s3} on books. How much money does she have left? Answer with only the final integer.", start - s1 - s2 - s3)
        # two-leg drive with pause trap
        speed, t1, speed2, t2, pause = rng.randint(38, 94), rng.randint(2, 6), rng.randint(32, 78), rng.randint(2, 5), rng.randint(1, 3)
        add("drive", f"A bus drives at {speed} km/h for {t1} hours, pauses {pause} hours at a depot, then drives {t2} more hours at {speed2} km/h. How many kilometers has it covered in total? Answer with only the final integer.", speed * t1 + speed2 * t2)
        # tray scaling with subtraction
        per, groups, extra, sold = rng.randint(9, 24), rng.randint(6, 17), rng.randint(11, 48), rng.randint(19, 88)
        add("bake", f"Each tray holds {per} rolls. A baker fills {groups} trays, adds {extra} loose rolls, and then sells {sold} rolls. How many rolls remain? Answer with only the final integer.", per * groups + extra - sold)
        # use/add/use ordering
        total, used, gift, used2 = rng.randint(280, 690), rng.randint(34, 96), rng.randint(21, 88), rng.randint(15, 74)
        add("tank", f"A tank holds {total} liters. Workers use {used} liters, a delivery adds {gift} liters, and then workers use {used2} more liters. How many liters are in the tank now? Answer with only the final integer.", total - used + gift - used2)
        # remainder inside a two-step
        stock, pack, sent = rng.randint(61, 190), rng.randint(6, 14), rng.randint(3, 9)
        add("packleft", f"A warehouse has {stock} lamps and packs them in crates of {pack}. It ships {sent} full crates. How many lamps remain unshipped? Answer with only the final integer.", stock - pack * sent)
        # four-actor relative chain
        x, y, z, w = rng.randint(24, 76), rng.randint(8, 34), rng.randint(6, 28), rng.randint(4, 22)
        add("cards", f"Ana has {x} cards. Ben has {y} more than Ana. Cara has {z} fewer than Ben. Dara has {w} more than Cara. How many cards does Dara have? Answer with only the final integer.", x + y - z + w)
    rng.shuffle(items)
    for i, it in enumerate(items):
        it["id"] = f"item-{i:04d}"
    return items


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20260930)
    ap.add_argument("--per-template", type=int, default=10)
    ap.add_argument("--out", default="items.json")
    a = ap.parse_args()
    its = gen_items(a.seed, a.per_template)
    with open(a.out, "w") as f:
        json.dump(its, f, indent=1)
    print(f"wrote {len(its)} items to {a.out}")
