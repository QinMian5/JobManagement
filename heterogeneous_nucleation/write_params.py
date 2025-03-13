# Author: Mian Qin
# Date Created: 9/23/24
import json


def main():
    x_star_list_list = [
        list(range(0, 130 + 1, 10)),
        list(range(140, 400 + 1, 20)),
        # [200, 220, 240, 260]
    ]
    kappa_list = [
        0.5,
        0.05,
    ]
    ramp_rate = 1000 / 3000  # 300K
    # ramp_rate = 1000 / 6000  # 270 K
    # ramp_rate = 1000 / 12000  # 250K
    prd_time = 10000
    job_params = {}
    x_star_init = 1090
    for x_star_list, kappa in zip(x_star_list_list, kappa_list):
        for x_star in x_star_list:
            delta_x_star = abs(x_star - x_star_init)
            ramp_time = int(delta_x_star / ramp_rate)
            nsteps = int((ramp_time + prd_time) / 0.002)
            job_params[f"op_{x_star}"] = {
                "QBAR": {"X_STAR": x_star, "X_STAR_INIT": x_star_init, "KAPPA": kappa, "PHI": 0.0},
                "TEMPERATURE": 300,
                "RAMP_TIME": ramp_time,
                "PRD_TIME": prd_time,
                "NSTEPS": nsteps
            }
    with open("job_params.json", 'w') as file:
        json.dump(job_params, file, indent='\t')


if __name__ == "__main__":
    main()
