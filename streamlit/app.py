import base64
import os
import sys
_dir = os.path.dirname(__file__)
sys.path.append(_dir)

import pandas as pd
import json
import re
import time
import uuid
from io import BytesIO
from pathlib import Path

import numpy as np
import streamlit as st
from PIL import Image
from streamlit_drawable_canvas import st_canvas
from svgpathtools import parse_path

import SessionState
import connector as cn


def main():
    handle = cn.Handler.getInstance()

    if 'button_id' not in st.session_state:
        st.session_state['button_id'] = ''
    if 'color_to_label' not in st.session_state:
        st.session_state['color_to_label'] = {}

    #color_annotation_app()
    new_calibration_app()

    with st.sidebar:
        st.markdown("---")
        st.markdown(
            '<h6>Made by Kelly @4D </h6>',
            unsafe_allow_html=True,
        )
        

def new_calibration_app():
    st.markdown(
        """
        MAKE LIFE A RIDE..
    """
    )
    groundtype = None
    obj = None
    handle = cn.Handler.getInstance()
    img_list = handle.bd.getImageList()
    bg_image = Image.open(img_list[0])
    gr_image = Image.open(handle.bd.gr_img)    
    bw, bh= bg_image.size[:2]    
    gw, gh = gr_image.size[:2]
    du_line_gr = []
    du_line_img = []

    print(bw, bh)
    print(gw, gh)

    label_color = (
        st.sidebar.color_picker("Annotation color: ", "#EA1010") + "77"
    )  # for alpha from 00 to FF
    label = st.sidebar.text_input("Label", "Default")

    canvas_result1 = st_canvas(
        fill_color=label_color,
        stroke_width=1,
        background_image=gr_image,
        drawing_mode='line',
        key="cvs1",
    )

    canvas_result2 = st_canvas(
        fill_color=label_color,
        stroke_width=1,
        background_image=bg_image,
        drawing_mode='line',
        key="cvs2",
    )
    if st.button('Select Complete.'):
        st.write('')
        handle.ExecuteExtract()
        time.sleep(10) 

    line_gr = []
    line_img = []
    if canvas_result1.json_data is not None :
        obj1 = canvas_result1.json_data['objects']
        for i in obj1 : 
            l = (i['left'] - i['width']/2) * (gw / 600)
            t = (i['top'] - i['height']/2) * (gh / 400)
            r = (i['left'] + i['width']/2) * (gw / 600)
            b = (i['top'] + i['height']/2) * (gh / 400)
            line_gr.append((l, t, r, b))

    if canvas_result2.json_data is not None :
        obj2 = canvas_result2.json_data['objects']
        for j in obj2 :
            l = (j['left'] - j['width']/2) * (bw / 600)
            t = (j['top'] - j['height']/2) * (bh / 400)
            r = (j['left'] + j['width']/2) * (bw / 600)
            b = (j['top'] + j['height']/2) * (bh / 400)

            line_img.append((l, t, r , b))

    du_line_gr = set(line_gr)
    du_line_img = set(line_img)
    handle.setRegion(du_line_gr, du_line_img)        


if __name__ == "__main__":
    st.set_page_config(
        page_title="EXODUSS", page_icon=":earthl2:"
    )
    st.title("CALIBRATION SIMULATION")
    st.sidebar.subheader("Config")
    main()
