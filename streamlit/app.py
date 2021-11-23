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
    playground = st.sidebar.selectbox("Type", options=list(handle.bd.getGroundType().keys()))
    
    print(playground)
    handle.setGround(playground)

    if 'button_id' not in st.session_state:
        st.session_state['button_id'] = ''
    if 'color_to_label' not in st.session_state:
        st.session_state['color_to_label'] = {}

    #color_annotation_app()
    center_circle_app()

    with st.sidebar:
        st.markdown("---")
        st.markdown(
            '<h6>Made by Kelly @4D </h6>',
            unsafe_allow_html=True,
        )
        

def color_annotation_app():
    st.markdown(
        """
        You can see one image.

        It is first Camera among camera group.
        
        Please Select ROI, 4 Region.
    """
    )
    groundtype = None
    scale = 5
    region = []
    obj = None
    handle = cn.Handler.getInstance()
    img_list = handle.bd.getImageList()
    bg_image = Image.open(img_list[0])
    w, h = bg_image.size[:2]
    bg_image = bg_image.resize((int(h/scale), int(w/scale)), Image.BILINEAR)
    print(w, h, bg_image.size[:2])

    label_color = (
        st.sidebar.color_picker("Annotation color: ", "#EA1010") + "77"
    )  # for alpha from 00 to FF
    label = st.sidebar.text_input("Label", "Default")

    canvas_result = st_canvas(
        fill_color=label_color,
        stroke_width=1,
        background_image=bg_image,
        height=h/scale,
        width=w/scale,
        drawing_mode='circle',
        key="color_annotation_app",
    )
    if st.button('Select Complete.'):
        st.write('')
        handle.ExecuteExtract()   
        #region.clear()     
        #canvas_result.json_data = None
        time.sleep(10) 

    if canvas_result.json_data is not None :
        obj = canvas_result.json_data['objects']

    if  obj is not None and len(obj)> 0:
        print(canvas_result.json_data['objects'])
        print("outer if .....")
        for ee in obj : 
            if ee['type'] == 'circle':
                df = pd.json_normalize(ee)
                region.append(df['top'], df['left'], df['width'])

        st.dataframe(df[["top", "left", "width", "height"]])        

            #print(region)
            #handle.setRegion(region)

#        st.session_state["color_to_label"][label_color] = label

        with st.expander("Color to label mapping"):
            st.json(st.session_state["color_to_label"])
    
    if st.button('Load Result'):
        path = "saved/"
        ilist = os.listdir(path)
        ilist.sort()
        for i in ilist:
            if 'save' not in i and 'png' in i:
                ii = Image.open(path + i)
                st.text(path + i)
                st.image(ii)


def center_circle_app():
    st.markdown(
        """
        You can see one image.
        It is first Camera among camera group.    
        Please select circle reiong for calibration.
    """
    )
    groundtype = None
    scale = 5
    region = []
    obj = None
    handle = cn.Handler.getInstance()
    img_list = handle.bd.getImageList()
    bg_image = Image.open(img_list[0])
    w, h = bg_image.size[:2]
    bg_image = bg_image.resize((int(h/scale), int(w/scale)), Image.BILINEAR)

    with open(_dir + "/saved_state.json", "r") as f:
        saved_state = json.load(f)

    canvas_result = st_canvas(
        fill_color="rgba(180, 165, 100, 0.3)",  # Fixed fill color with some opacity
        stroke_width=2,
        stroke_color="black",
        background_image=bg_image,
        initial_drawing=saved_state
        if st.sidebar.checkbox("Initialize with saved state", False)
        else None,
        height=h/scale,
        width=w/scale,
        drawing_mode="circle",
        key="center_circle_app",
    )
    #with st.echo("below"):
    if canvas_result.json_data is not None:
        df = pd.json_normalize(canvas_result.json_data["objects"])
        if len(df) == 0:
            return
        df["center_x"] = df["left"] + df["radius"] * np.cos(
            df["angle"] * np.pi / 180
        )
        df["center_y"] = df["top"] + df["radius"] * np.sin(
            df["angle"] * np.pi / 180
        )

        st.subheader("List of circle drawings")
        for _, row in df.iterrows():
            cx = row["center_x"] * scale
            cy = row["center_y"] * scale
            r = row["radius"] * scale
            st.markdown(
                f'Center coords: ({cx:.2f}, {cy:.2f}). Radius: {r:.2f}'
            )
            region.append(cx)
            region.append(cy)
            region.append(r)


    if st.button('Select Complete.'):
        handle.setRegion(region)        
        handle.ExecuteExtract()  


if __name__ == "__main__":
    st.set_page_config(
        page_title="GENESIS", page_icon=":earthl2:"
    )
    st.title("CALIBRATION SIMULATION")
    st.sidebar.subheader("Config")
    main()
