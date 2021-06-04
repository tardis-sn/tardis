.. raw:: html

    <div class = "twitter-block"  style = "width:100%; display: flex; justify-content: center;">
    <a class="twitter-timeline" data-width="850" data-height="900" href="https://twitter.com/tardis_sn?ref_src=twsrc%5Etfw">Tweets by tardis_sn</a> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>
    </div>



    <script>
        document.getElementsByClassName("twitter-block")[0].addEventListener(
          "DOMSubtreeModified",
          function (e) {
            for (var target = e.target; target && target != this; target = target.parentNode) {
              // loop parent nodes from the target to the delegation node
              if (target.matches("#twitter-widget-0")) {
                customizeTweetMedia();
              }
            }
          },
          false
        );

        var customizeTweetMedia = function () {
          var iframe = document.getElementById("twitter-widget-0");

          var cssLink = iframe.contentWindow.document.createElement("link");
          cssLink.setAttribute("rel", "stylesheet");
          cssLink.setAttribute("type", "text/css");
          cssLink.setAttribute("href", "./_static/css/theme.css");

          iframeHead = iframe.contentWindow.document.getElementsByTagName("head")[0];

          if (!iframeHead.contains(iframe.contentDocument.getElementById("sphinx-css"))) {
            iframeHead.appendChild(cssLink);
          }

          iframe.contentWindow.document.body.style.overflow = "hidden";

          var headerTitle = iframe.contentWindow.document.getElementsByClassName("timeline-Header-title")[0];
          headerTitle.style.fontSize = "18px";

          var headerTitleByInline = iframe.contentWindow.document.getElementsByClassName("timeline-Header-byline")[0];
          headerTitleByInline.style.fontSize = "18.5px";

          var tweetTextList = iframe.contentWindow.document.getElementsByClassName("timeline-Tweet-text");
          var index;
          for (let index = 0; index < tweetTextList.length; index++) {
            tweetTextList[index].style.fontSize = "16px";
            tweetTextList[index].style.fontFamily = "Lato, proxima-nova, Helvetica Neue, Arial, sans-serif";
            tweetTextList[index].style.lineHeight = "1.1";
          }

          tweetList = iframe.contentWindow.document.getElementsByClassName("timeline-TweetList")[0];
          tweetList.addEventListener(
            "DOMSubtreeModified",
            function () {
              customizeTweetMedia();
            },
            false
          );
        };

        // checks if the widget hasn't loaded and prints a message after 4 seconds
        setTimeout(() => {
          if (document.querySelector("#twitter-widget-0") == null) {
            const fallbackDiv = document.createElement("div");
            const linkSpan = document.createElement("span");
            const link = document.createElement("a");

            link.setAttribute("href", `https://twitter.com/tardis_sn`);
            link.textContent = "TARDIS Tweets";

            linkSpan.appendChild(link);

            const messageP = document.createElement("P");
            messageP.textContent =
              "It seems Twitter Widget was not loaded. Please check your internet connection or consider disabling tracking protection if on Firefox. If the problem persists, please contact us.";

            fallbackDiv.appendChild(linkSpan);
            fallbackDiv.appendChild(messageP);

            fallbackDiv.style.cssText = "font-family: Lato, sans-serif; text-align:center;";
            document.querySelector(".twitter-timeline").replaceWith(fallbackDiv);
          }
        }, 4000);

    </script>