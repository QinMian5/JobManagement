# Author: Mian Qin
# Date Created: 9/23/24
import json


def main():
    icing = True
    x_star_list = list(range(100, 1701, 100))
    # ramp_rate = 200 / 10000
    prd_time = 3000
    job_params = {}
    x_star_init = 0 if icing else 2000
    for x_star in x_star_list:
        delta_x_star = abs(x_star - x_star_init)
        # ramp_time = int(delta_x_star / ramp_rate)
        ramp_time = 2000
        nsteps = int((ramp_time + prd_time) / 0.002)
        job_params[f"op_{x_star}"] = {
            "QBAR": {"TYPE": "parabola", "STAR": x_star, "STAR_INIT": x_star_init, "KAPPA": 0.0029},
            "TEMPERATURE": 300,
            "RAMP_TIME": ramp_time,
            "PRD_TIME": prd_time,
            "NSTEPS": nsteps
        }
    with open("job_params.json", 'w') as file:
        json.dump(job_params, file, indent='\t')


if __name__ == "__main__":
    main()
