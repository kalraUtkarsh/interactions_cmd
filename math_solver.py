import selenium
from selenium import webdriver
from selenium.webdriver.common.by import By 

wd = webdriver.Chrome('./chromedriver.exe')

initial_url = """https://tbot.xyz/math/#eyJ1Ijo1NjMwOTM3MTEsIm4iOiJVdGthcnNoIEthbHJhIiwiZyI6Ik1hdGhCYXR0bGUiLCJjaSI6Ii0xNjY3MDExMjY0NzA4MTIyMzAyIiwiaSI6IkJRQUFBQWowQWdEUElKQWhxQjRMUGpNUjdaSSJ9ZjI0OGFkZjAyNDE5MmE3ZmNkZmRiZGNmOWM4OTFmN2Y=&tgShareScoreUrl=tg%3A%2F%2Fshare_game_score%3Fhash%3Df2xxcp5wFqHB-AAkVDCZPMA_DdhvaXWgaRvy0aC0asB9Se0mSBnV8ezkg6YiCGZv"""
wd.get(initial_url)

wd.find_element(By.ID, "button_correct").click()
for b in range(57):
    x = wd.find_element(By.ID,"task_x").text
    y = wd.find_element(By.ID,"task_y").text
    op = wd.find_element(By.ID,"task_op").text
    res = wd.find_element(By.ID,"task_res").text

    equation = f"{x} {op} {y}"
    equation = equation.replace('×', '*').replace('–','-')
    real_res = int(eval(equation))

    if real_res == int(res):
        wd.find_element(By.ID, "button_correct").click()
    else:
        wd.find_element(By.ID, "button_wrong").click()