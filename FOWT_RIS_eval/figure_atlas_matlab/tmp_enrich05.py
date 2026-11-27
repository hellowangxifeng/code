p = "fig05_link.m"
s = open(p, encoding="utf-8").read()
blocks = [
    ("S.pt_tx", "pt_OP"),
    ("S.pt_tx", "pt_ABER"),
    ("S.dist_R", "dist_OP"),
    ("S.inc_theta", "inc_OP"),
    ("S.ap_radius", "ap_OP"),
    ("S.wind_U", "wind_OP"),
    ("S.winddir_x", "winddir_OP"),
    ("S.winddir_x", "winddir_ABER"),
]
for x, pre in blocks:
    old = (f"semilogy_smooth(ax,{x},S.{pre}_linear_c151,C_LIN,'-');\n"
           f"semilogy_smooth(ax,{x},S.{pre}_focus_c151,C_FOC,'-');")
    new = old + (f"\nsemilogy_smooth(ax,{x},S.{pre}_linear_c25,C_LIN,'-.');"
                 f"\nsemilogy_smooth(ax,{x},S.{pre}_focus_c25,C_FOC,'-.');")
    if old not in s:
        print("MISS", pre)
    s = s.replace(old, new)
open(p, "w", encoding="utf-8").write(s)
print("done")
