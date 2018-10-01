import cPickle, base64

try:
    from SimpleSession.versions.v65 import (
        beginRestore,
        registerAfterModelsCB,
        reportRestoreError,
        checkVersion,
    )
except ImportError:
    from chimera import UserError

    raise UserError(
        "Cannot open session that was saved in a"
        " newer version of Chimera; update your version"
    )
checkVersion([1, 12, 41613])
import chimera
from chimera import replyobj

replyobj.status("Restoring session...", blankAfter=0)
replyobj.status("Beginning session restore...", blankAfter=0, secondary=True)
beginRestore()


def restoreCoreModels():
    from SimpleSession.versions.v65 import (
        init,
        restoreViewer,
        restoreMolecules,
        restoreColors,
        restoreSurfaces,
        restoreVRML,
        restorePseudoBondGroups,
        restoreModelAssociations,
    )

    molInfo = cPickle.loads(
        base64.b64decode(
            "gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwFOfYdVCWJhbGxTY2FsZXEDSwFHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLAUc/8AAAAAAAAH2HVQVjb2xvcnEFSwFLAH2HVQpyaWJib25UeXBlcQZLAUsAfYdVCnN0aWNrU2NhbGVxB0sBRz/wAAAAAAAAfYdVDG1tQ0lGSGVhZGVyc3EIXXEJTmFVDGFyb21hdGljTW9kZXEKSwFLAX2HVQp2ZHdEZW5zaXR5cQtLAUdAFAAAAAAAAH2HVQZoaWRkZW5xDEsBiX2HVQ1hcm9tYXRpY0NvbG9ycQ1LAU59h1UPcmliYm9uU21vb3RoaW5ncQ5LAUsAfYdVCWF1dG9jaGFpbnEPSwGIfYdVCnBkYlZlcnNpb25xEEsBSwB9h1UIb3B0aW9uYWxxEX1VD2xvd2VyQ2FzZUNoYWluc3ESSwGJfYdVCWxpbmVXaWR0aHETSwFHP/AAAAAAAAB9h1UPcmVzaWR1ZUxhYmVsUG9zcRRLAUsAfYdVBG5hbWVxFUsBWBcAAABhbGlnbmVkX3dpdGhfZHVtbXkucnN0N32HVQ9hcm9tYXRpY0Rpc3BsYXlxFksBiX2HVQ9yaWJib25TdGlmZm5lc3NxF0sBRz/pmZmZmZmafYdVCnBkYkhlYWRlcnNxGF1xGX1xGmFVA2lkc3EbSwFLAEsAhn2HVQ5zdXJmYWNlT3BhY2l0eXEcSwFHv/AAAAAAAAB9h1UQYXJvbWF0aWNMaW5lVHlwZXEdSwFLAn2HVRRyaWJib25IaWRlc01haW5jaGFpbnEeSwGIfYdVB2Rpc3BsYXlxH0sBiH2HdS4="
        )
    )
    resInfo = cPickle.loads(
        base64.b64decode(
            "gAJ9cQEoVQZpbnNlcnRxAksFVQEgfYdVC2ZpbGxEaXNwbGF5cQNLBYl9h1UEbmFtZXEESwVYAwAAAERNMX1xBShYAwAAAERNMl1xBksDYVgDAAAARE0zXXEHSwRhWAMAAABDQjZdcQhLAGFYAwAAAEJVVF1xCUsBYXWHVQVjaGFpbnEKSwVYAQAAACB9h1UOcmliYm9uRHJhd01vZGVxC0sFSwJ9h1UCc3NxDEsFiYmGfYdVCG1vbGVjdWxlcQ1LBUsAfYdVC3JpYmJvbkNvbG9ycQ5LBUsDfXEPKEsBTl1xEEsASwGGcRFhhksCTl1xEksBSwGGcRNhhnWHVQVsYWJlbHEUSwVYAAAAAH2HVQpsYWJlbENvbG9ycRVLBUsDfXEWKEsBTl1xF0sASwGGcRhhhksCTl1xGUsBSwGGcRphhnWHVQhmaWxsTW9kZXEbSwVLAX2HVQVpc0hldHEcSwWJfYdVC2xhYmVsT2Zmc2V0cR1LBU59h1UIcG9zaXRpb25xHl1xH0sBSwWGcSBhVQ1yaWJib25EaXNwbGF5cSFLBYl9h1UIb3B0aW9uYWxxIn1VBHNzSWRxI0sFSv////99h3Uu"
        )
    )
    atomInfo = cPickle.loads(
        base64.b64decode(
            "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"
        )
    )
    bondInfo = cPickle.loads(
        base64.b64decode(
            "gAJ9cQEoVQVjb2xvcnECS4tOfYdVBWF0b21zcQNdcQQoXXEFKEtrS21lXXEGKEtqS2xlXXEHKEtnS2llXXEIKEtnS2plXXEJKEtmS2hlXXEKKEtmS2tlXXELKEtjS2plXXEMKEtiS2tlXXENKEthS2ZlXXEOKEthS2dlXXEPKEtgS2FlXXEQKEtgS2JlXXERKEtgS2NlXXESKEtZS1tlXXETKEtYS1plXXEUKEtVS1dlXXEVKEtVS1hlXXEWKEtUS1ZlXXEXKEtUS1llXXEYKEtRS1hlXXEZKEtRS2llXXEaKEtQS1llXXEbKEtQS2hlXXEcKEtPS1RlXXEdKEtPS1VlXXEeKEtOS09lXXEfKEtOS1BlXXEgKEtOS1FlXXEhKEtHS0llXXEiKEtGS0hlXXEjKEtDS0VlXXEkKEtDS0ZlXXElKEtCS0RlXXEmKEtCS0dlXXEnKEs/S0ZlXXEoKEs/S1dlXXEpKEs+S0dlXXEqKEs+S1ZlXXErKEs9S0JlXXEsKEs9S0NlXXEtKEs8Sz1lXXEuKEs8Sz5lXXEvKEs8Sz9lXXEwKEs1SzdlXXExKEs0SzZlXXEyKEsxSzNlXXEzKEsxSzRlXXE0KEswSzJlXXE1KEswSzVlXXE2KEstSzRlXXE3KEstS0VlXXE4KEssSzVlXXE5KEssS0RlXXE6KEsrSzBlXXE7KEsrSzFlXXE8KEsqSytlXXE9KEsqSyxlXXE+KEsqSy1lXXE/KEsjSyVlXXFAKEsiSyRlXXFBKEsfSyFlXXFCKEsfSyJlXXFDKEseSyBlXXFEKEseSyNlXXFFKEsbSyJlXXFGKEsbSzNlXXFHKEsaSyNlXXFIKEsaSzJlXXFJKEsZSx5lXXFKKEsZSx9lXXFLKEsYSxllXXFMKEsYSxplXXFNKEsYSxtlXXFOKEsRSxNlXXFPKEsQSxJlXXFQKEsPS2NlXXFRKEsOS2JlXXFSKEsNSw9lXXFTKEsNSxBlXXFUKEsMSw5lXXFVKEsMSxFlXXFWKEsJSxBlXXFXKEsJSyFlXXFYKEsISxFlXXFZKEsISyBlXXFaKEsHSwxlXXFbKEsHSw1lXXFcKEsGSwdlXXFdKEsGSwhlXXFeKEsGSwllXXFfKEt4S3xlXXFgKEt2S3hlXXFhKEtyS3ZlXXFiKEtpS3BlXXFjKEtpS3FlXXFkKEtoS25lXXFlKEtoS29lXXFmKEthS2VlXXFnKEtgS2RlXXFoKEtXS15lXXFpKEtXS19lXXFqKEtWS1xlXXFrKEtWS11lXXFsKEtPS1NlXXFtKEtOS1JlXXFuKEtFS0xlXXFvKEtFS01lXXFwKEtES0plXXFxKEtES0tlXXFyKEs9S0FlXXFzKEs8S0BlXXF0KEszSzplXXF1KEszSztlXXF2KEsySzhlXXF3KEsySzllXXF4KEsrSy9lXXF5KEsqSy5lXXF6KEshSyhlXXF7KEshSyllXXF8KEsgSyZlXXF9KEsgSydlXXF+KEsZSx1lXXF/KEsYSxxlXXGAKEsPSxZlXXGBKEsPSxdlXXGCKEsOSxRlXXGDKEsOSxVlXXGEKEsHSwtlXXGFKEsGSwplXXGGKEt8S31lXXGHKEt8S35lXXGIKEt8S39lXXGJKEt4S3plXXGKKEt4S3tlXXGLKEt2S3dlXXGMKEt2S3llXXGNKEtyS3NlXXGOKEtyS3RlXXGPKEtyS3VlZVUFbGFiZWxxkEuLWAAAAAB9h1UIaGFsZmJvbmRxkUuLiH2HVQZyYWRpdXNxkkuLRz/JmZmgAAAAfYdVC2xhYmVsT2Zmc2V0cZNLi059h1UIZHJhd01vZGVxlEuLSwF9h1UIb3B0aW9uYWxxlX1VB2Rpc3BsYXlxlkuLSwJ9h3Uu"
        )
    )
    crdInfo = cPickle.loads(
        base64.b64decode(
            "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"
        )
    )
    surfInfo = {
        "category": (0, None, {}),
        "probeRadius": (0, None, {}),
        "pointSize": (0, None, {}),
        "name": [],
        "density": (0, None, {}),
        "colorMode": (0, None, {}),
        "useLighting": (0, None, {}),
        "transparencyBlendMode": (0, None, {}),
        "molecule": [],
        "smoothLines": (0, None, {}),
        "lineWidth": (0, None, {}),
        "allComponents": (0, None, {}),
        "twoSidedLighting": (0, None, {}),
        "customVisibility": [],
        "drawMode": (0, None, {}),
        "display": (0, None, {}),
        "customColors": [],
    }
    vrmlInfo = {
        "subid": (0, None, {}),
        "display": (0, None, {}),
        "id": (0, None, {}),
        "vrmlString": [],
        "name": (0, None, {}),
    }
    colors = {
        u"Ru": ((0.141176, 0.560784, 0.560784), 1, u"default"),
        u"Re": ((0.14902, 0.490196, 0.670588), 1, u"default"),
        u"Rf": ((0.8, 0, 0.34902), 1, u"default"),
        u"Ra": ((0, 0.490196, 0), 1, u"default"),
        u"Rb": ((0.439216, 0.180392, 0.690196), 1, u"default"),
        u"Rn": ((0.258824, 0.509804, 0.588235), 1, u"default"),
        u"Rh": ((0.0392157, 0.490196, 0.54902), 1, u"default"),
        u"Be": ((0.760784, 1, 0), 1, u"default"),
        u"Ba": ((0, 0.788235, 0), 1, u"default"),
        u"Bh": ((0.878431, 0, 0.219608), 1, u"default"),
        u"Bi": ((0.619608, 0.309804, 0.709804), 1, u"default"),
        u"Bk": ((0.541176, 0.309804, 0.890196), 1, u"default"),
        u"Br": ((0.65098, 0.160784, 0.160784), 1, u"default"),
        u"H": ((1, 1, 1), 1, u"default"),
        u"P": ((1, 0.501961, 0), 1, u"default"),
        u"Os": ((0.14902, 0.4, 0.588235), 1, u"default"),
        u"Ge": ((0.4, 0.560784, 0.560784), 1, u"default"),
        u"Gd": ((0.270588, 1, 0.780392), 1, u"default"),
        u"Ga": ((0.760784, 0.560784, 0.560784), 1, u"default"),
        u"Pr": ((0.85098, 1, 0.780392), 1, u"default"),
        u"Pt": ((0.815686, 0.815686, 0.878431), 1, u"default"),
        u"Pu": ((0, 0.419608, 1), 1, u"default"),
        u"C": ((0.564706, 0.564706, 0.564706), 1, u"default"),
        u"grey": ((0.745098, 0.745098, 0.745098), 1, u"default"),
        u"Pb": ((0.341176, 0.34902, 0.380392), 1, u"default"),
        u"Pa": ((0, 0.631373, 1), 1, u"default"),
        u"Pd": ((0, 0.411765, 0.521569), 1, u"default"),
        u"Xe": ((0.258824, 0.619608, 0.690196), 1, u"default"),
        u"Po": ((0.670588, 0.360784, 0), 1, u"default"),
        u"Pm": ((0.639216, 1, 0.780392), 1, u"default"),
        u"Hs": ((0.901961, 0, 0.180392), 1, u"default"),
        u"Ho": ((0, 1, 0.611765), 1, u"default"),
        u"Hf": ((0.301961, 0.760784, 1), 1, u"default"),
        u"Hg": ((0.721569, 0.721569, 0.815686), 1, u"default"),
        u"He": ((0.85098, 1, 1), 1, u"default"),
        u"Md": ((0.701961, 0.0509804, 0.65098), 1, u"default"),
        u"Mg": ((0.541176, 1, 0), 1, u"default"),
        u"K": ((0.560784, 0.25098, 0.831373), 1, u"default"),
        u"Mn": ((0.611765, 0.478431, 0.780392), 1, u"default"),
        u"O": ((1, 0.0509804, 0.0509804), 1, u"default"),
        u"Zr": ((0.580392, 0.878431, 0.878431), 1, u"default"),
        u"S": ((1, 1, 0.188235), 1, u"default"),
        u"W": ((0.129412, 0.580392, 0.839216), 1, u"default"),
        u"Zn": ((0.490196, 0.501961, 0.690196), 1, u"default"),
        u"Mt": ((0.921569, 0, 0.14902), 1, u"default"),
        u"Eu": ((0.380392, 1, 0.780392), 1, u"default"),
        u"Es": ((0.701961, 0.121569, 0.831373), 1, u"default"),
        u"Er": ((0, 0.901961, 0.458824), 1, u"default"),
        u"Ni": ((0.313725, 0.815686, 0.313725), 1, u"default"),
        u"No": ((0.741176, 0.0509804, 0.529412), 1, u"default"),
        u"Na": ((0.670588, 0.360784, 0.94902), 1, u"default"),
        u"Nb": ((0.45098, 0.760784, 0.788235), 1, u"default"),
        u"Nd": ((0.780392, 1, 0.780392), 1, u"default"),
        u"Ne": ((0.701961, 0.890196, 0.960784), 1, u"default"),
        u"Np": ((0, 0.501961, 1), 1, u"default"),
        u"Fr": ((0.258824, 0, 0.4), 1, u"default"),
        u"Fe": ((0.878431, 0.4, 0.2), 1, u"default"),
        u"Fm": ((0.701961, 0.121569, 0.729412), 1, u"default"),
        u"B": ((1, 0.709804, 0.709804), 1, u"default"),
        u"F": ((0.564706, 0.878431, 0.313725), 1, u"default"),
        u"Sr": ((0, 1, 0), 1, u"default"),
        u"N": ((0.188235, 0.313725, 0.972549), 1, u"default"),
        u"Kr": ((0.360784, 0.721569, 0.819608), 1, u"default"),
        u"Si": ((0.941176, 0.784314, 0.627451), 1, u"default"),
        u"Sn": ((0.4, 0.501961, 0.501961), 1, u"default"),
        u"Sm": ((0.560784, 1, 0.780392), 1, u"default"),
        u"V": ((0.65098, 0.65098, 0.670588), 1, u"default"),
        u"Sc": ((0.901961, 0.901961, 0.901961), 1, u"default"),
        u"Sb": ((0.619608, 0.388235, 0.709804), 1, u"default"),
        u"Sg": ((0.85098, 0, 0.270588), 1, u"default"),
        u"Se": ((1, 0.631373, 0), 1, u"default"),
        u"Co": ((0.941176, 0.564706, 0.627451), 1, u"default"),
        u"Cm": ((0.470588, 0.360784, 0.890196), 1, u"default"),
        u"Cl": ((0.121569, 0.941176, 0.121569), 1, u"default"),
        u"Ca": ((0.239216, 1, 0), 1, u"default"),
        u"Cf": ((0.631373, 0.211765, 0.831373), 1, u"default"),
        u"Ce": ((1, 1, 0.780392), 1, u"default"),
        u"Cd": ((1, 0.85098, 0.560784), 1, u"default"),
        u"Lu": ((0, 0.670588, 0.141176), 1, u"default"),
        u"Cs": ((0.341176, 0.0901961, 0.560784), 1, u"default"),
        u"Cr": ((0.541176, 0.6, 0.780392), 1, u"default"),
        u"Cu": ((0.784314, 0.501961, 0.2), 1, u"default"),
        u"La": ((0.439216, 0.831373, 1), 1, u"default"),
        u"Li": ((0.8, 0.501961, 1), 1, u"default"),
        u"Tl": ((0.65098, 0.329412, 0.301961), 1, u"default"),
        u"Tm": ((0, 0.831373, 0.321569), 1, u"default"),
        u"Lr": ((0.780392, 0, 0.4), 1, u"default"),
        u"Th": ((0, 0.729412, 1), 1, u"default"),
        u"Ti": ((0.74902, 0.760784, 0.780392), 1, u"default"),
        u"Te": ((0.831373, 0.478431, 0), 1, u"default"),
        u"Tb": ((0.188235, 1, 0.780392), 1, u"default"),
        u"Tc": ((0.231373, 0.619608, 0.619608), 1, u"default"),
        u"Ta": ((0.301961, 0.65098, 1), 1, u"default"),
        u"Yb": ((0, 0.74902, 0.219608), 1, u"default"),
        u"Db": ((0.819608, 0, 0.309804), 1, u"default"),
        u"Dy": ((0.121569, 1, 0.780392), 1, u"default"),
        u"I": ((0.580392, 0, 0.580392), 1, u"default"),
        u"U": ((0, 0.560784, 1), 1, u"default"),
        u"Y": ((0.580392, 1, 1), 1, u"default"),
        u"Ac": ((0.439216, 0.670588, 0.980392), 1, u"default"),
        u"Ag": ((0.752941, 0.752941, 0.752941), 1, u"default"),
        u"Ir": ((0.0901961, 0.329412, 0.529412), 1, u"default"),
        u"Am": ((0.329412, 0.360784, 0.94902), 1, u"default"),
        u"Al": ((0.74902, 0.65098, 0.65098), 1, u"default"),
        u"As": ((0.741176, 0.501961, 0.890196), 1, u"default"),
        u"Ar": ((0.501961, 0.819608, 0.890196), 1, u"default"),
        u"Au": ((1, 0.819608, 0.137255), 1, u"default"),
        u"At": ((0.458824, 0.309804, 0.270588), 1, u"default"),
        u"In": ((0.65098, 0.458824, 0.45098), 1, u"default"),
        u"Mo": ((0.329412, 0.709804, 0.709804), 1, u"default"),
    }
    materials = {u"default": ((0, 0, 0), 30)}
    pbInfo = {
        "category": [u"distance monitor"],
        "bondInfo": [
            {
                "color": (0, None, {}),
                "atoms": [],
                "label": (0, None, {}),
                "halfbond": (0, None, {}),
                "labelColor": (0, None, {}),
                "labelOffset": (0, None, {}),
                "drawMode": (0, None, {}),
                "display": (0, None, {}),
            }
        ],
        "lineType": (1, 2, {}),
        "color": (1, 130, {}),
        "optional": {"fixedLabels": (True, False, (1, False, {}))},
        "display": (1, True, {}),
        "showStubBonds": (1, False, {}),
        "lineWidth": (1, 1, {}),
        "stickScale": (1, 1, {}),
        "id": [-2],
    }
    modelAssociations = {}
    colorInfo = (
        133,
        (u"", (0.745098, 0.745098, 0.745098, 0.2)),
        {
            (u"green", (0, 1, 0, 1)): [2],
            (u"", (0, 0, 1, 0.2)): [4, 58, 95],
            (u"", (1, 0, 0, 0.2)): [126, 127, 128],
            (u"", (0.106829, 0.702586, 0.652042, 1)): [0],
            (u"", (1, 1, 1, 1)): [131],
            (u"grey", (0.745098, 0.745098, 0.745098, 1)): [129],
            (u"red", (1, 0, 0, 1)): [3],
            (u"", (0, 1, 0, 0.2)): [112, 122],
            (u"", (0.545455, 0, 1, 1)): [132],
            (u"blue", (0, 0, 1, 1)): [1],
            (u"yellow", (1, 1, 0, 1)): [130],
        },
    )
    viewerInfo = {
        "cameraAttrs": {
            "center": (0.031, -0.0945, -3.0070000190735),
            "fieldOfView": 17.183565421784,
            "nearFar": (7.2154204221993, -13.229420460346),
            "ortho": True,
            "eyeSeparation": 50.8,
            "focal": -3.0070000190735,
        },
        "viewerAttrs": {
            "silhouetteColor": None,
            "clipping": False,
            "showSilhouette": True,
            "showShadows": False,
            "viewSize": 11.391911504425,
            "labelsOnTop": True,
            "depthCueRange": (0.5, 1),
            "silhouetteWidth": 2,
            "singleLayerTransparency": True,
            "shadowTextureSize": 2048,
            "backgroundImage": [None, 1, 2, 1, 0, 0],
            "backgroundGradient": [
                ("Chimera default", [(1, 1, 1, 1), (0, 0, 1, 1)], 1),
                1,
                0,
                0,
            ],
            "depthCue": True,
            "highlight": 0,
            "scaleFactor": 1.1936205300669,
            "angleDependentTransparency": True,
            "backgroundMethod": 0,
        },
        "viewerHL": 132,
        "cameraMode": "mono",
        "detail": 1.5,
        "viewerFog": None,
        "viewerBG": 131,
    }

    replyobj.status("Initializing session restore...", blankAfter=0, secondary=True)
    from SimpleSession.versions.v65 import expandSummary

    init(dict(enumerate(expandSummary(colorInfo))))
    replyobj.status("Restoring colors...", blankAfter=0, secondary=True)
    restoreColors(colors, materials)
    replyobj.status("Restoring molecules...", blankAfter=0, secondary=True)
    restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
    replyobj.status("Restoring surfaces...", blankAfter=0, secondary=True)
    restoreSurfaces(surfInfo)
    replyobj.status("Restoring VRML models...", blankAfter=0, secondary=True)
    restoreVRML(vrmlInfo)
    replyobj.status("Restoring pseudobond groups...", blankAfter=0, secondary=True)
    restorePseudoBondGroups(pbInfo)
    replyobj.status("Restoring model associations...", blankAfter=0, secondary=True)
    restoreModelAssociations(modelAssociations)
    replyobj.status("Restoring camera...", blankAfter=0, secondary=True)
    restoreViewer(viewerInfo)


try:
    restoreCoreModels()
except:
    reportRestoreError("Error restoring core models")

    replyobj.status("Restoring extension info...", blankAfter=0, secondary=True)


try:
    import StructMeasure
    from StructMeasure.DistMonitor import restoreDistances

    registerAfterModelsCB(restoreDistances, 1)
except:
    reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
    formattedPositions = {}
    import Midas

    Midas.restoreMidasBase(formattedPositions)


try:
    restoreMidasBase()
except:
    reportRestoreError("Error restoring Midas base state")


def restoreMidasText():
    from Midas import midas_text

    midas_text.aliases = {}
    midas_text.userSurfCategories = {}


try:
    restoreMidasText()
except:
    reportRestoreError("Error restoring Midas text state")


def restore_cap_attributes():
    cap_attributes = {
        "cap_attributes": [],
        "cap_color": None,
        "cap_offset": 0.01,
        "class": "Caps_State",
        "default_cap_offset": 0.01,
        "mesh_style": False,
        "shown": True,
        "subdivision_factor": 1.0,
        "version": 1,
    }
    import SurfaceCap.session

    SurfaceCap.session.restore_cap_attributes(cap_attributes)


registerAfterModelsCB(restore_cap_attributes)


def restore_volume_data():
    volume_data_state = {
        "class": "Volume_Manager_State",
        "data_and_regions_state": [],
        "version": 2,
    }
    from VolumeViewer import session

    session.restore_volume_data_state(volume_data_state)


try:
    restore_volume_data()
except:
    reportRestoreError("Error restoring volume data")

geomData = {"AxisManager": {}, "CentroidManager": {}, "PlaneManager": {}}

try:
    from StructMeasure.Geometry import geomManager

    geomManager._restoreSession(geomData)
except:
    reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
    import SimpleSession
    import RibbonStyleEditor

    userScalings = []
    userXSections = []
    userResidueClasses = []
    residueData = [
        (1, "Chimera default", "rounded", u"unknown"),
        (2, "Chimera default", "rounded", u"unknown"),
        (3, "Chimera default", "rounded", u"unknown"),
        (4, "Chimera default", "rounded", u"unknown"),
        (5, "Chimera default", "rounded", u"unknown"),
    ]
    flags = RibbonStyleEditor.NucleicDefault1
    SimpleSession.registerAfterModelsCB(
        RibbonStyleEditor.restoreState,
        (userScalings, userXSections, userResidueClasses, residueData, flags),
    )


try:
    restoreSession_RibbonStyleEditor()
except:
    reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = "gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4="
scPickle = "gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=="
kfPickle = (
    "gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu"
)


def restoreAnimation():
    "A method to unpickle and restore animation objects"
    # Scenes must be unpickled after restoring transitions, because each
    # scene links to a 'scene' transition. Likewise, keyframes must be
    # unpickled after restoring scenes, because each keyframe links to a scene.
    # The unpickle process is left to the restore* functions, it's
    # important that it doesn't happen prior to calling those functions.
    import SimpleSession
    from Animate.Session import restoreTransitions
    from Animate.Session import restoreScenes
    from Animate.Session import restoreKeyframes

    SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
    SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
    SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)


try:
    restoreAnimation()
except:
    reportRestoreError("Error in Animate.Session")


def restoreLightController():
    import Lighting

    Lighting._setFromParams(
        {
            "ratio": 1.25,
            "brightness": 1.16,
            "material": [30.0, (0.85, 0.85, 0.85), 0.0],
            "back": [
                (0.3574067443365933, 0.6604015517481455, -0.6604015517481456),
                (1.0, 1.0, 1.0),
                0.0,
            ],
            "mode": "two-point",
            "key": [
                (-0.3574067443365933, 0.6604015517481455, 0.6604015517481456),
                (1.0, 1.0, 1.0),
                1.0,
            ],
            "contrast": 0.83,
            "fill": [
                (0.2505628070857316, 0.2505628070857316, 0.9351131265310294),
                (1.0, 1.0, 1.0),
                0.0,
            ],
        }
    )


try:
    restoreLightController()
except:
    reportRestoreError("Error restoring lighting parameters")

mdmovieData = {
    "length": 1,
    "startFrame": None,
    "endFrame": None,
    "molecule": 0,
    "name": "aligned_with_dummy.rst7",
}

try:
    from Movie import restoreSession

    mdmovie = restoreSession(mdmovieData)
except:
    reportRestoreError("Error restoring MD Movie interface")


def restoreRemainder():
    from SimpleSession.versions.v65 import (
        restoreWindowSize,
        restoreOpenStates,
        restoreSelections,
        restoreFontInfo,
        restoreOpenModelsAttrs,
        restoreModelClip,
        restoreSilhouettes,
    )

    curSelIds = []
    savedSels = []
    openModelsAttrs = {"cofrMethod": 4}
    windowSize = (888, 678)
    xformMap = {
        0: (
            (
                (-0.035831583629827, 0.99807456861732, -0.050628579813384),
                89.214908418703,
            ),
            (3.0382954392039, -0.039825362408825, -2.9354095863633),
            True,
        )
    }
    fontInfo = {"face": ("Sans Serif", "Normal", 16)}
    clipPlaneInfo = {}
    silhouettes = {0: True, 270: True}

    replyobj.status("Restoring window...", blankAfter=0, secondary=True)
    restoreWindowSize(windowSize)
    replyobj.status("Restoring open states...", blankAfter=0, secondary=True)
    restoreOpenStates(xformMap)
    replyobj.status("Restoring font info...", blankAfter=0, secondary=True)
    restoreFontInfo(fontInfo)
    replyobj.status("Restoring selections...", blankAfter=0, secondary=True)
    restoreSelections(curSelIds, savedSels)
    replyobj.status("Restoring openModel attributes...", blankAfter=0, secondary=True)
    restoreOpenModelsAttrs(openModelsAttrs)
    replyobj.status("Restoring model clipping...", blankAfter=0, secondary=True)
    restoreModelClip(clipPlaneInfo)
    replyobj.status("Restoring per-model silhouettes...", blankAfter=0, secondary=True)
    restoreSilhouettes(silhouettes)

    replyobj.status(
        "Restoring remaining extension info...", blankAfter=0, secondary=True
    )


try:
    restoreRemainder()
except:
    reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs

makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore

replyobj.status("Finishing restore...", blankAfter=0, secondary=True)
endRestore({})
replyobj.status("", secondary=True)
replyobj.status("Restore finished.")
