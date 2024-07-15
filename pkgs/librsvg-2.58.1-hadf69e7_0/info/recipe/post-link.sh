GDK_PIXBUF_POST_LINK_SCRIPT="${PREFIX}/bin/.gdk-pixbuf-post-link.sh"
if [ -x "$GDK_PIXBUF_POST_LINK_SCRIPT" ]
then
    # The gdk-pixbuf post-link function updates the loaders
    "$GDK_PIXBUF_POST_LINK_SCRIPT"
else
    cat >> ${PREFIX}/.messages.txt << EOF
librsvg: The post-link script did not complete.
To take advantage of gdk-pixbuf's support for librsvg, please run:
    $GDK_PIXBUF_POST_LINK_SCRIPT
EOF
fi
